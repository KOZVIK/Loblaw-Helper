"""Initialize the SQLite database and load cell-count.csv.

Run from the repository root (or any other working directory) with:

    python load_data.py

Only Python's standard library is required.
"""

from __future__ import annotations

import csv
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Final


ROOT: Final = Path(__file__).resolve().parent
CSV_PATH: Final = ROOT / "cell-count.csv"
DB_PATH: Final = ROOT / "cell-count.db"

POPULATIONS: Final = (
    "b_cell",
    "cd8_t_cell",
    "cd4_t_cell",
    "nk_cell",
    "monocyte",
)

REQUIRED_COLUMNS: Final = {
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
}

SCHEMA: Final = """
PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS cell_populations;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id   INTEGER PRIMARY KEY,
    project_name TEXT NOT NULL UNIQUE
);

CREATE TABLE subjects (
    subject_id   INTEGER PRIMARY KEY,
    project_id   INTEGER NOT NULL,
    subject_code TEXT NOT NULL,
    condition    TEXT NOT NULL,
    age          INTEGER NOT NULL CHECK (age >= 0),
    gender       TEXT NOT NULL,
    treatment    TEXT NOT NULL,
    response     TEXT CHECK (response IN ('yes', 'no') OR response IS NULL),
    FOREIGN KEY (project_id) REFERENCES projects(project_id),
    UNIQUE (project_id, subject_code)
);

CREATE TABLE samples (
    sample_id                  INTEGER PRIMARY KEY,
    subject_id                 INTEGER NOT NULL,
    sample_code                TEXT NOT NULL UNIQUE,
    sample_type                TEXT NOT NULL,
    time_from_treatment_start  INTEGER NOT NULL,
    FOREIGN KEY (subject_id) REFERENCES subjects(subject_id)
);

CREATE TABLE cell_populations (
    population_id   INTEGER PRIMARY KEY,
    population_name TEXT NOT NULL UNIQUE
);

CREATE TABLE cell_counts (
    sample_id     INTEGER NOT NULL,
    population_id INTEGER NOT NULL,
    cell_count    INTEGER NOT NULL CHECK (cell_count >= 0),
    PRIMARY KEY (sample_id, population_id),
    FOREIGN KEY (sample_id) REFERENCES samples(sample_id),
    FOREIGN KEY (population_id) REFERENCES cell_populations(population_id)
);

CREATE INDEX idx_subjects_condition_treatment_response
    ON subjects(condition, treatment, response);
CREATE INDEX idx_samples_subject_time_type
    ON samples(subject_id, time_from_treatment_start, sample_type);
CREATE INDEX idx_cell_counts_population
    ON cell_counts(population_id);
"""


def required_int(row: dict[str, str], column: str, line_number: int) -> int:
    """Parse a required integer and report the source line on failure."""
    value = row[column].strip()
    try:
        return int(value)
    except ValueError as exc:
        raise ValueError(
            f"CSV line {line_number}: {column!r} must be an integer, got {value!r}"
        ) from exc


def load_database(
    csv_path: Path = CSV_PATH, db_path: Path = DB_PATH
) -> tuple[int, int]:
    """Rebuild the database from the CSV and return subject/sample counts."""
    if not csv_path.is_file():
        raise FileNotFoundError(f"Input file not found: {csv_path}")

    with csv_path.open(newline="", encoding="utf-8-sig") as source:
        reader = csv.DictReader(source)
        columns = set(reader.fieldnames or ())
        missing = REQUIRED_COLUMNS - columns
        if missing:
            raise ValueError(
                "Input CSV is missing required columns: " + ", ".join(sorted(missing))
            )
        rows = list(reader)

    with closing(sqlite3.connect(db_path)) as connection:
        with connection:
            connection.executescript(SCHEMA)
            connection.executemany(
                "INSERT INTO cell_populations (population_name) VALUES (?)",
                ((name,) for name in POPULATIONS),
            )
            population_ids = dict(
                connection.execute(
                    "SELECT population_name, population_id FROM cell_populations"
                )
            )

            project_ids: dict[str, int] = {}
            subject_ids: dict[tuple[str, str], int] = {}
            subject_metadata: dict[tuple[str, str], tuple[object, ...]] = {}

            for line_number, row in enumerate(rows, start=2):
                project_name = row["project"].strip()
                subject_code = row["subject"].strip()
                subject_key = (project_name, subject_code)

                if project_name not in project_ids:
                    cursor = connection.execute(
                        "INSERT INTO projects (project_name) VALUES (?)",
                        (project_name,),
                    )
                    project_ids[project_name] = cursor.lastrowid

                response = row["response"].strip().lower() or None
                metadata = (
                    row["condition"].strip(),
                    required_int(row, "age", line_number),
                    row["sex"].strip(),
                    row["treatment"].strip(),
                    response,
                )

                if subject_key not in subject_ids:
                    cursor = connection.execute(
                        """
                        INSERT INTO subjects (
                            project_id, subject_code, condition, age, gender,
                            treatment, response
                        ) VALUES (?, ?, ?, ?, ?, ?, ?)
                        """,
                        (project_ids[project_name], subject_code, *metadata),
                    )
                    subject_ids[subject_key] = cursor.lastrowid
                    subject_metadata[subject_key] = metadata
                elif subject_metadata[subject_key] != metadata:
                    raise ValueError(
                        f"CSV line {line_number}: inconsistent metadata for "
                        f"subject {subject_code!r} in project {project_name!r}"
                    )

                sample_cursor = connection.execute(
                    """
                    INSERT INTO samples (
                        subject_id, sample_code, sample_type,
                        time_from_treatment_start
                    ) VALUES (?, ?, ?, ?)
                    """,
                    (
                        subject_ids[subject_key],
                        row["sample"].strip(),
                        row["sample_type"].strip(),
                        required_int(row, "time_from_treatment_start", line_number),
                    ),
                )
                sample_id = sample_cursor.lastrowid

                connection.executemany(
                    """
                    INSERT INTO cell_counts (sample_id, population_id, cell_count)
                    VALUES (?, ?, ?)
                    """,
                    (
                        (
                            sample_id,
                            population_ids[population],
                            required_int(row, population, line_number),
                        )
                        for population in POPULATIONS
                    ),
                )

            subject_count = connection.execute(
                "SELECT COUNT(*) FROM subjects"
            ).fetchone()[0]
            sample_count = connection.execute(
                "SELECT COUNT(*) FROM samples"
            ).fetchone()[0]

    return subject_count, sample_count


def main() -> None:
    subject_count, sample_count = load_database()
    print(
        f"Created {DB_PATH.name} with "
        f"{subject_count:,} subjects and {sample_count:,} samples."
    )


if __name__ == "__main__":
    main()
