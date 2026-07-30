"""Part 4: query baseline clinical-trial subsets and summary counts."""

from __future__ import annotations

import csv
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Final


ROOT: Final = Path(__file__).resolve().parent
DB_PATH: Final = ROOT / "cell-count.db"
OUTPUT_DIR: Final = ROOT / "outputs"
BASELINE_SAMPLES_PATH: Final = OUTPUT_DIR / "baseline_melanoma_miraclib_pbmc.csv"
SUBSET_COUNTS_PATH: Final = OUTPUT_DIR / "baseline_subset_counts.csv"
B_CELL_ANSWER_PATH: Final = OUTPUT_DIR / "baseline_male_responder_b_cell_average.txt"

BASELINE_FILTER: Final = """
subjects.condition = 'melanoma'
AND subjects.treatment = 'miraclib'
AND samples.sample_type = 'PBMC'
AND samples.time_from_treatment_start = 0
"""

BASELINE_SAMPLES_QUERY: Final = f"""
SELECT
    projects.project_name AS project,
    subjects.subject_code AS subject,
    subjects.condition,
    subjects.age,
    subjects.gender,
    subjects.treatment,
    subjects.response,
    samples.sample_code AS sample,
    samples.sample_type,
    samples.time_from_treatment_start
FROM samples
JOIN subjects USING (subject_id)
JOIN projects USING (project_id)
WHERE {BASELINE_FILTER}
ORDER BY projects.project_name, subjects.subject_code, samples.sample_code
"""

PROJECT_COUNTS_QUERY: Final = f"""
SELECT projects.project_name, COUNT(*) AS sample_count
FROM samples
JOIN subjects USING (subject_id)
JOIN projects USING (project_id)
WHERE {BASELINE_FILTER}
GROUP BY projects.project_id, projects.project_name
ORDER BY projects.project_name
"""

RESPONSE_COUNTS_QUERY: Final = f"""
SELECT subjects.response, COUNT(DISTINCT subjects.subject_id) AS subject_count
FROM samples
JOIN subjects USING (subject_id)
WHERE {BASELINE_FILTER}
GROUP BY subjects.response
ORDER BY subjects.response
"""

GENDER_COUNTS_QUERY: Final = f"""
SELECT subjects.gender, COUNT(DISTINCT subjects.subject_id) AS subject_count
FROM samples
JOIN subjects USING (subject_id)
WHERE {BASELINE_FILTER}
GROUP BY subjects.gender
ORDER BY subjects.gender
"""

B_CELL_AVERAGE_QUERY: Final = """
SELECT AVG(cell_counts.cell_count)
FROM samples
JOIN subjects USING (subject_id)
JOIN cell_counts USING (sample_id)
JOIN cell_populations USING (population_id)
WHERE subjects.condition = 'melanoma'
  AND subjects.gender = 'M'
  AND subjects.response = 'yes'
  AND samples.time_from_treatment_start = 0
  AND cell_populations.population_name = 'b_cell'
"""

BASELINE_COLUMNS: Final = (
    "project",
    "subject",
    "condition",
    "age",
    "gender",
    "treatment",
    "response",
    "sample",
    "sample_type",
    "time_from_treatment_start",
)


def query_part4(db_path: Path = DB_PATH) -> dict[str, object]:
    """Return all Part 4 query results."""
    if not db_path.is_file():
        raise FileNotFoundError(
            f"Database not found: {db_path}. Run `python load_data.py` first."
        )

    with closing(sqlite3.connect(db_path)) as connection:
        baseline_samples = connection.execute(BASELINE_SAMPLES_QUERY).fetchall()
        project_counts = dict(connection.execute(PROJECT_COUNTS_QUERY))
        response_counts = dict(connection.execute(RESPONSE_COUNTS_QUERY))
        gender_counts = dict(connection.execute(GENDER_COUNTS_QUERY))
        b_cell_average = connection.execute(B_CELL_AVERAGE_QUERY).fetchone()[0]

    if b_cell_average is None:
        raise ValueError("No samples matched the B-cell average query")

    return {
        "baseline_samples": baseline_samples,
        "project_counts": project_counts,
        "response_counts": response_counts,
        "gender_counts": gender_counts,
        "b_cell_average": float(b_cell_average),
    }


def write_outputs(
    results: dict[str, object],
    baseline_path: Path = BASELINE_SAMPLES_PATH,
    counts_path: Path = SUBSET_COUNTS_PATH,
    answer_path: Path = B_CELL_ANSWER_PATH,
) -> None:
    """Write the sample list, grouped counts, and requested average."""
    baseline_path.parent.mkdir(parents=True, exist_ok=True)

    with baseline_path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(BASELINE_COLUMNS)
        writer.writerows(results["baseline_samples"])

    with counts_path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(("metric", "category", "count"))
        for project, count in results["project_counts"].items():
            writer.writerow(("samples_by_project", project, count))
        for response, count in results["response_counts"].items():
            writer.writerow(("subjects_by_response", response, count))
        for gender, count in results["gender_counts"].items():
            writer.writerow(("subjects_by_gender", gender, count))

    answer_path.write_text(
        "Average B-cell count for responding melanoma males at baseline "
        f"(all sample and treatment types): {results['b_cell_average']:.2f}\n",
        encoding="utf-8",
    )


def main() -> None:
    results = query_part4()
    write_outputs(results)

    print(
        "Baseline melanoma/miraclib/PBMC samples: "
        f"{len(results['baseline_samples']):,}"
    )
    print("\nSamples by project:")
    for project, count in results["project_counts"].items():
        print(f"  {project}: {count}")
    print("\nDistinct subjects by response:")
    for response, count in results["response_counts"].items():
        print(f"  {response}: {count}")
    print("\nDistinct subjects by gender:")
    for gender, count in results["gender_counts"].items():
        print(f"  {gender}: {count}")
    print(
        "\nAverage B-cell count for responding melanoma males at baseline "
        f"(all sample and treatment types): {results['b_cell_average']:.2f}"
    )
    print(f"\nWrote sample list to {BASELINE_SAMPLES_PATH.relative_to(ROOT)}.")
    print(f"Wrote subset counts to {SUBSET_COUNTS_PATH.relative_to(ROOT)}.")
    print(f"Wrote B-cell answer to {B_CELL_ANSWER_PATH.relative_to(ROOT)}.")


if __name__ == "__main__":
    main()
