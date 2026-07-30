"""Tests for the Part 1 CSV-to-SQLite loader."""

from __future__ import annotations

import csv
import sqlite3
import tempfile
import unittest
from pathlib import Path

from analyze_frequencies import OUTPUT_COLUMNS, frequency_rows, write_summary
from load_data import POPULATIONS, load_database


FIELDNAMES = [
    "project",
    "subject",
    "condition",
    "age",
    "sex",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
    *POPULATIONS,
]


def sample_row(**overrides: str) -> dict[str, str]:
    row = {
        "project": "prj1",
        "subject": "sbj001",
        "condition": "melanoma",
        "age": "57",
        "sex": "M",
        "treatment": "miraclib",
        "response": "yes",
        "sample": "sample001",
        "sample_type": "PBMC",
        "time_from_treatment_start": "0",
        "b_cell": "10",
        "cd8_t_cell": "20",
        "cd4_t_cell": "30",
        "nk_cell": "40",
        "monocyte": "50",
    }
    row.update(overrides)
    return row


class LoadDataTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_directory.cleanup)
        self.root = Path(self.temp_directory.name)
        self.csv_path = self.root / "cell-count.csv"
        self.db_path = self.root / "cell-count.db"

    def write_csv(
        self,
        rows: list[dict[str, str]],
        fieldnames: list[str] = FIELDNAMES,
    ) -> None:
        with self.csv_path.open("w", newline="", encoding="utf-8") as output:
            writer = csv.DictWriter(output, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)

    def connect(self) -> sqlite3.Connection:
        connection = sqlite3.connect(self.db_path)
        self.addCleanup(connection.close)
        return connection

    def test_loads_normalized_schema_and_data(self) -> None:
        self.write_csv(
            [
                sample_row(),
                sample_row(
                    sample="sample002",
                    time_from_treatment_start="7",
                    b_cell="12",
                ),
            ]
        )

        counts = load_database(self.csv_path, self.db_path)
        connection = self.connect()

        self.assertEqual(counts, (1, 2))
        table_counts = {
            table: connection.execute(
                f"SELECT COUNT(*) FROM {table}"
            ).fetchone()[0]
            for table in (
                "projects",
                "subjects",
                "samples",
                "cell_populations",
                "cell_counts",
            )
        }
        self.assertEqual(
            table_counts,
            {
                "projects": 1,
                "subjects": 1,
                "samples": 2,
                "cell_populations": 5,
                "cell_counts": 10,
            },
        )
        subject = connection.execute(
            "SELECT condition, gender, treatment, response FROM subjects"
        ).fetchone()
        self.assertEqual(subject, ("melanoma", "M", "miraclib", "yes"))
        self.assertEqual(
            connection.execute(
                """
                SELECT cell_count
                FROM cell_counts
                JOIN samples USING (sample_id)
                JOIN cell_populations USING (population_id)
                WHERE sample_code = 'sample002'
                  AND population_name = 'b_cell'
                """
            ).fetchone()[0],
            12,
        )
        self.assertEqual(connection.execute("PRAGMA foreign_key_check").fetchall(), [])

    def test_blank_response_is_stored_as_null(self) -> None:
        self.write_csv([sample_row(response="")])

        load_database(self.csv_path, self.db_path)
        response = self.connect().execute(
            "SELECT response FROM subjects"
        ).fetchone()[0]

        self.assertIsNone(response)

    def test_rerun_rebuilds_database_without_duplicate_rows(self) -> None:
        self.write_csv([sample_row()])

        load_database(self.csv_path, self.db_path)
        load_database(self.csv_path, self.db_path)
        connection = self.connect()

        self.assertEqual(
            connection.execute("SELECT COUNT(*) FROM samples").fetchone()[0], 1
        )
        self.assertEqual(
            connection.execute("SELECT COUNT(*) FROM cell_counts").fetchone()[0], 5
        )

    def test_rejects_missing_required_column(self) -> None:
        fieldnames = [name for name in FIELDNAMES if name != "condition"]
        row = sample_row()
        row.pop("condition")
        self.write_csv([row], fieldnames)

        with self.assertRaisesRegex(ValueError, "condition"):
            load_database(self.csv_path, self.db_path)

    def test_rejects_inconsistent_subject_metadata(self) -> None:
        self.write_csv(
            [
                sample_row(),
                sample_row(sample="sample002", age="58"),
            ]
        )

        with self.assertRaisesRegex(ValueError, "inconsistent metadata"):
            load_database(self.csv_path, self.db_path)

    def test_rejects_invalid_integer(self) -> None:
        self.write_csv([sample_row(b_cell="not-a-number")])

        with self.assertRaisesRegex(ValueError, "'b_cell' must be an integer"):
            load_database(self.csv_path, self.db_path)


class FrequencyAnalysisTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_directory.cleanup)
        self.root = Path(self.temp_directory.name)
        self.csv_path = self.root / "cell-count.csv"
        self.db_path = self.root / "cell-count.db"
        self.output_path = self.root / "outputs" / "cell_frequencies.csv"

    def write_csv(self, rows: list[dict[str, str]]) -> None:
        with self.csv_path.open("w", newline="", encoding="utf-8") as output:
            writer = csv.DictWriter(output, fieldnames=FIELDNAMES)
            writer.writeheader()
            writer.writerows(rows)

    def test_calculates_total_and_population_percentages(self) -> None:
        self.write_csv([sample_row()])
        load_database(self.csv_path, self.db_path)

        rows = frequency_rows(self.db_path)

        self.assertEqual(len(rows), 5)
        self.assertEqual(
            rows[0],
            ("sample001", 150, "b_cell", 10, 100.0 * 10 / 150),
        )
        self.assertAlmostEqual(sum(float(row[4]) for row in rows), 100.0)

    def test_zero_total_produces_zero_percentages(self) -> None:
        zero_counts = {population: "0" for population in POPULATIONS}
        self.write_csv([sample_row(**zero_counts)])
        load_database(self.csv_path, self.db_path)

        rows = frequency_rows(self.db_path)

        self.assertTrue(all(row[1] == 0 for row in rows))
        self.assertTrue(all(row[4] == 0.0 for row in rows))

    def test_writes_required_csv_columns_and_values(self) -> None:
        self.write_csv([sample_row()])
        load_database(self.csv_path, self.db_path)

        count = write_summary(frequency_rows(self.db_path), self.output_path)

        with self.output_path.open(newline="", encoding="utf-8") as source:
            output_rows = list(csv.DictReader(source))
        self.assertEqual(count, 5)
        self.assertEqual(list(output_rows[0]), list(OUTPUT_COLUMNS))
        self.assertEqual(output_rows[0]["sample"], "sample001")
        self.assertEqual(output_rows[0]["total_count"], "150")
        self.assertEqual(output_rows[0]["population"], "b_cell")
        self.assertEqual(output_rows[0]["count"], "10")
        self.assertEqual(output_rows[0]["percentage"], "6.666667")


if __name__ == "__main__":
    unittest.main()
