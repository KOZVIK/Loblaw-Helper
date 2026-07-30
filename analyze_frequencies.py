"""Generate the Part 2 cell-population relative-frequency summary."""

from __future__ import annotations

import csv
import sqlite3
from contextlib import closing
from pathlib import Path
from typing import Final, Iterable


ROOT: Final = Path(__file__).resolve().parent
DB_PATH: Final = ROOT / "cell-count.db"
OUTPUT_PATH: Final = ROOT / "outputs" / "cell_frequencies.csv"
OUTPUT_COLUMNS: Final = (
    "sample",
    "total_count",
    "population",
    "count",
    "percentage",
)

FREQUENCY_QUERY: Final = """
WITH sample_totals AS (
    SELECT sample_id, SUM(cell_count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
)
SELECT
    samples.sample_code AS sample,
    sample_totals.total_count,
    cell_populations.population_name AS population,
    cell_counts.cell_count AS count,
    CASE
        WHEN sample_totals.total_count = 0 THEN 0.0
        ELSE 100.0 * cell_counts.cell_count / sample_totals.total_count
    END AS percentage
FROM samples
JOIN sample_totals USING (sample_id)
JOIN cell_counts USING (sample_id)
JOIN cell_populations USING (population_id)
ORDER BY samples.sample_code, cell_populations.population_id
"""


def frequency_rows(db_path: Path = DB_PATH) -> list[tuple[object, ...]]:
    """Return one relative-frequency row per sample and cell population."""
    if not db_path.is_file():
        raise FileNotFoundError(
            f"Database not found: {db_path}. Run `python load_data.py` first."
        )

    with closing(sqlite3.connect(db_path)) as connection:
        return connection.execute(FREQUENCY_QUERY).fetchall()


def write_summary(
    rows: Iterable[tuple[object, ...]], output_path: Path = OUTPUT_PATH
) -> int:
    """Write frequency rows to CSV and return the number written."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    row_count = 0
    with output_path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output)
        writer.writerow(OUTPUT_COLUMNS)
        for sample, total_count, population, count, percentage in rows:
            writer.writerow(
                (
                    sample,
                    total_count,
                    population,
                    count,
                    f"{float(percentage):.6f}",
                )
            )
            row_count += 1
    return row_count


def print_preview(rows: list[tuple[object, ...]], limit: int = 10) -> None:
    """Print a compact preview without requiring third-party packages."""
    formatted_rows = [
        (*row[:4], f"{float(row[4]):.6f}") for row in rows[:limit]
    ]
    preview = [OUTPUT_COLUMNS, *formatted_rows]
    widths = [
        max(len(str(row[index])) for row in preview)
        for index in range(len(OUTPUT_COLUMNS))
    ]
    for row_number, row in enumerate(preview):
        print(
            "  ".join(
                str(value).ljust(widths[index])
                for index, value in enumerate(row)
            )
        )
        if row_number == 0:
            print("  ".join("-" * width for width in widths))


def main() -> None:
    rows = frequency_rows()
    row_count = write_summary(rows)
    print_preview(rows)
    print(f"\nWrote {row_count:,} rows to {OUTPUT_PATH.relative_to(ROOT)}.")


if __name__ == "__main__":
    main()
