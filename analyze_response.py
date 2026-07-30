"""Part 3: compare immune-cell frequencies by miraclib response."""

from __future__ import annotations

import csv
import os
import sqlite3
from contextlib import closing
from pathlib import Path
from statistics import median
from typing import Final, Sequence

ROOT: Final = Path(__file__).resolve().parent
os.environ.setdefault("MPLCONFIGDIR", str(ROOT / ".matplotlib-cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

from load_data import POPULATIONS


DB_PATH: Final = ROOT / "cell-count.db"
OUTPUT_DIR: Final = ROOT / "outputs"
STATISTICS_PATH: Final = OUTPUT_DIR / "response_statistics.csv"
REPORT_PATH: Final = OUTPUT_DIR / "response_statistics.txt"
BOXPLOT_PATH: Final = OUTPUT_DIR / "response_boxplots.png"
ALPHA: Final = 0.05

SUBJECT_FREQUENCY_QUERY: Final = """
WITH sample_totals AS (
    SELECT sample_id, SUM(cell_count) AS total_count
    FROM cell_counts
    GROUP BY sample_id
),
sample_frequencies AS (
    SELECT
        subjects.subject_id,
        subjects.response,
        cell_populations.population_name AS population,
        CASE
            WHEN sample_totals.total_count = 0 THEN 0.0
            ELSE 100.0 * cell_counts.cell_count / sample_totals.total_count
        END AS percentage
    FROM subjects
    JOIN samples USING (subject_id)
    JOIN sample_totals USING (sample_id)
    JOIN cell_counts USING (sample_id)
    JOIN cell_populations USING (population_id)
    WHERE subjects.condition = 'melanoma'
      AND subjects.treatment = 'miraclib'
      AND samples.sample_type = 'PBMC'
      AND subjects.response IN ('yes', 'no')
)
SELECT
    subject_id,
    response,
    population,
    AVG(percentage) AS mean_percentage
FROM sample_frequencies
GROUP BY subject_id, response, population
ORDER BY population, response, subject_id
"""

STATISTICS_COLUMNS: Final = (
    "population",
    "responder_subjects",
    "non_responder_subjects",
    "responder_median_percentage",
    "non_responder_median_percentage",
    "median_difference",
    "mann_whitney_u",
    "p_value",
    "adjusted_p_value",
    "significant",
    "rank_biserial_correlation",
)


def benjamini_hochberg(p_values: Sequence[float]) -> list[float]:
    """Return Benjamini–Hochberg false-discovery-rate adjusted p-values."""
    count = len(p_values)
    if count == 0:
        return []
    order = sorted(range(count), key=p_values.__getitem__)
    adjusted = [0.0] * count
    running_minimum = 1.0
    for rank_index in range(count - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        candidate = min(1.0, p_values[original_index] * count / rank)
        running_minimum = min(running_minimum, candidate)
        adjusted[original_index] = running_minimum
    return adjusted


def load_subject_frequencies(
    db_path: Path = DB_PATH,
) -> dict[str, dict[str, list[float]]]:
    """Load subject-level mean percentages for each response and population."""
    if not db_path.is_file():
        raise FileNotFoundError(
            f"Database not found: {db_path}. Run `python load_data.py` first."
        )

    groups = {
        population: {"yes": [], "no": []} for population in POPULATIONS
    }
    with closing(sqlite3.connect(db_path)) as connection:
        rows = connection.execute(SUBJECT_FREQUENCY_QUERY).fetchall()
    for _subject_id, response, population, percentage in rows:
        groups[population][response].append(float(percentage))
    return groups


def calculate_statistics(
    groups: dict[str, dict[str, list[float]]],
) -> list[dict[str, object]]:
    """Calculate two-sided tests, effect sizes, and FDR-adjusted p-values."""
    results: list[dict[str, object]] = []
    raw_p_values: list[float] = []

    for population in POPULATIONS:
        responders = groups[population]["yes"]
        non_responders = groups[population]["no"]
        if not responders or not non_responders:
            raise ValueError(
                f"Both response groups are required for population {population!r}"
            )

        test = mannwhitneyu(
            responders, non_responders, alternative="two-sided", method="auto"
        )
        responder_median = median(responders)
        non_responder_median = median(non_responders)
        rank_biserial = (
            2.0 * float(test.statistic) / (len(responders) * len(non_responders))
            - 1.0
        )
        raw_p_values.append(float(test.pvalue))
        results.append(
            {
                "population": population,
                "responder_subjects": len(responders),
                "non_responder_subjects": len(non_responders),
                "responder_median_percentage": responder_median,
                "non_responder_median_percentage": non_responder_median,
                "median_difference": responder_median - non_responder_median,
                "mann_whitney_u": float(test.statistic),
                "p_value": float(test.pvalue),
                "adjusted_p_value": 0.0,
                "significant": False,
                "rank_biserial_correlation": rank_biserial,
            }
        )

    for result, adjusted_p_value in zip(
        results, benjamini_hochberg(raw_p_values)
    ):
        result["adjusted_p_value"] = adjusted_p_value
        result["significant"] = adjusted_p_value < ALPHA
    return results


def write_statistics(
    results: list[dict[str, object]],
    statistics_path: Path = STATISTICS_PATH,
    report_path: Path = REPORT_PATH,
) -> None:
    """Write machine-readable statistics and a concise text conclusion."""
    statistics_path.parent.mkdir(parents=True, exist_ok=True)
    with statistics_path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=STATISTICS_COLUMNS)
        writer.writeheader()
        for result in results:
            formatted = dict(result)
            for column in (
                "responder_median_percentage",
                "non_responder_median_percentage",
                "median_difference",
                "mann_whitney_u",
                "p_value",
                "adjusted_p_value",
                "rank_biserial_correlation",
            ):
                formatted[column] = f"{float(result[column]):.6g}"
            formatted["significant"] = (
                "yes" if result["significant"] else "no"
            )
            writer.writerow(formatted)

    significant = [
        str(result["population"]) for result in results if result["significant"]
    ]
    conclusion = (
        "Significant populations after FDR correction: "
        + (", ".join(significant) if significant else "none")
        + "."
    )
    report_path.write_text(
        "\n".join(
            (
                "Miraclib response analysis: melanoma PBMC samples",
                "",
                "Statistical unit: each subject's mean population frequency "
                "across available PBMC timepoints.",
                "Test: two-sided Mann-Whitney U for responders versus "
                "non-responders.",
                "Multiple testing: Benjamini-Hochberg FDR across five "
                f"populations; significance threshold q < {ALPHA}.",
                conclusion,
                "",
            )
        ),
        encoding="utf-8",
    )


def create_boxplots(
    groups: dict[str, dict[str, list[float]]],
    output_path: Path = BOXPLOT_PATH,
) -> None:
    """Create one responder/non-responder boxplot per immune population."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, len(POPULATIONS), figsize=(16, 5), sharey=True)
    colors = ("#2878B5", "#D64A3A")

    for axis, population in zip(axes, POPULATIONS):
        values = [groups[population]["yes"], groups[population]["no"]]
        plot = axis.boxplot(
            values,
            tick_labels=["Responders", "Non-\nresponders"],
            patch_artist=True,
            widths=0.6,
            medianprops={"color": "black", "linewidth": 1.5},
        )
        for patch, color in zip(plot["boxes"], colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.75)
        axis.set_title(population.replace("_", " ").title())
        axis.grid(axis="y", alpha=0.25)

    axes[0].set_ylabel("Subject mean relative frequency (%)")
    figure.suptitle(
        "Melanoma PBMC cell frequencies by miraclib response",
        fontsize=14,
        fontweight="bold",
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    groups = load_subject_frequencies()
    results = calculate_statistics(groups)
    write_statistics(results)
    create_boxplots(groups)

    print("population       responders  non-responders  p-value     FDR q-value  significant")
    print("----------------  ----------  --------------  ----------  -----------  -----------")
    for result in results:
        print(
            f"{str(result['population']):16}"
            f"  {int(result['responder_subjects']):10}"
            f"  {int(result['non_responder_subjects']):14}"
            f"  {float(result['p_value']):10.4g}"
            f"  {float(result['adjusted_p_value']):11.4g}"
            f"  {'yes' if result['significant'] else 'no':>11}"
        )
    print(f"\nWrote statistics to {STATISTICS_PATH.relative_to(ROOT)}.")
    print(f"Wrote conclusions to {REPORT_PATH.relative_to(ROOT)}.")
    print(f"Wrote boxplots to {BOXPLOT_PATH.relative_to(ROOT)}.")


if __name__ == "__main__":
    main()
