"""Tests for the Part 4 subset queries."""

from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from analyze_subsets import query_part4, write_outputs
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


def row(**overrides: str) -> dict[str, str]:
    result = {
        "project": "prj1",
        "subject": "sbj001",
        "condition": "melanoma",
        "age": "50",
        "sex": "M",
        "treatment": "miraclib",
        "response": "yes",
        "sample": "sample001",
        "sample_type": "PBMC",
        "time_from_treatment_start": "0",
        "b_cell": "100",
        "cd8_t_cell": "20",
        "cd4_t_cell": "30",
        "nk_cell": "40",
        "monocyte": "50",
    }
    result.update(overrides)
    return result


class SubsetAnalysisTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp_directory.cleanup)
        self.root = Path(self.temp_directory.name)
        self.csv_path = self.root / "cell-count.csv"
        self.db_path = self.root / "cell-count.db"

    def build_database(self) -> None:
        rows = [
            row(),
            row(
                project="prj2",
                subject="sbj002",
                sample="sample002",
                response="no",
                sex="F",
                b_cell="200",
            ),
            row(
                project="prj1",
                subject="sbj003",
                sample="sample003",
                treatment="phauximab",
                b_cell="300",
            ),
            row(
                project="prj1",
                subject="sbj004",
                sample="sample004",
                condition="carcinoma",
                b_cell="900",
            ),
        ]
        with self.csv_path.open("w", newline="", encoding="utf-8") as output:
            writer = csv.DictWriter(output, fieldnames=FIELDNAMES)
            writer.writeheader()
            writer.writerows(rows)
        load_database(self.csv_path, self.db_path)

    def test_baseline_subset_and_grouped_counts(self) -> None:
        self.build_database()

        results = query_part4(self.db_path)

        self.assertEqual(len(results["baseline_samples"]), 2)
        self.assertEqual(results["project_counts"], {"prj1": 1, "prj2": 1})
        self.assertEqual(results["response_counts"], {"no": 1, "yes": 1})
        self.assertEqual(results["gender_counts"], {"F": 1, "M": 1})

    def test_b_cell_average_includes_all_treatments(self) -> None:
        self.build_database()

        results = query_part4(self.db_path)

        self.assertEqual(results["b_cell_average"], 200.0)

    def test_writes_part4_outputs(self) -> None:
        self.build_database()
        results = query_part4(self.db_path)
        baseline_path = self.root / "outputs" / "baseline.csv"
        counts_path = self.root / "outputs" / "counts.csv"
        answer_path = self.root / "outputs" / "answer.txt"

        write_outputs(results, baseline_path, counts_path, answer_path)

        with baseline_path.open(newline="", encoding="utf-8") as source:
            baseline_rows = list(csv.DictReader(source))
        with counts_path.open(newline="", encoding="utf-8") as source:
            count_rows = list(csv.reader(source))
        self.assertEqual(len(baseline_rows), 2)
        self.assertIn(
            ["subjects_by_gender", "M", "1"],
            count_rows,
        )
        self.assertIn("200.00", answer_path.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
