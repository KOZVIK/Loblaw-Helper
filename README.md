# Loblaw Bio Immune Cell Analysis

Python and SQLite project for analyzing immune cell populations in a clinical
trial.

## Quick start

In GitHub Codespaces, run these commands from the repository root:

```bash
make setup
make pipeline
make dashboard
```

- `make setup` installs all Python dependencies.
- `make pipeline` rebuilds the database and runs Parts 1–4.
- `make dashboard` starts the interactive Streamlit dashboard.

Open the forwarded Streamlit port when Codespaces prompts you.

## Outputs

The pipeline creates `cell-count.db` and the following files in `outputs/`:

| File | Contents |
| --- | --- |
| `cell_frequencies.csv` | Cell count and relative frequency for every sample and population |
| `response_boxplots.png` | Miraclib responder versus non-responder boxplots |
| `response_statistics.csv` | Part 3 statistical results |
| `response_statistics.txt` | Short Part 3 conclusion |
| `baseline_melanoma_miraclib_pbmc.csv` | Part 4 baseline sample subset |
| `baseline_subset_counts.csv` | Counts by project, response, and gender |
| `baseline_male_responder_b_cell_average.txt` | Requested B-cell average |

## Main results

- The database contains 3 projects, 3,500 subjects, 10,500 samples, and
  52,500 population-count records.
- Part 3 compares 331 responders with 325 non-responders. It uses each
  subject's mean PBMC frequency, two-sided Mann–Whitney U tests, and
  Benjamini–Hochberg correction across five populations.
- No population is significant at an adjusted p-value below 0.05. CD4 T cells
  have the smallest raw p-value (`p = 0.0124`, adjusted `q = 0.0621`).
- Part 4 identifies 656 baseline melanoma PBMC samples treated with miraclib:
  384 from `prj1` and 272 from `prj3`.
- The baseline subset contains 331 responder and 325 non-responder subjects;
  312 are female and 344 are male.
- The average B-cell count for responding melanoma males at time 0, across all
  sample and treatment types, is **10206.15**.

## Database schema

```text
projects 1 --< subjects 1 --< samples 1 --< cell_counts >-- 1 cell_populations
```

| Table | Purpose |
| --- | --- |
| `projects` | One row per project |
| `subjects` | Subject condition, age, gender, treatment, and response |
| `samples` | Sample ID, type, and time from treatment start |
| `cell_populations` | Population-name lookup table |
| `cell_counts` | Count for each sample and population |

This normalized design avoids repeating project, subject, and sample metadata
for every cell population. Foreign keys protect relationships, uniqueness
constraints prevent duplicate identifiers, and checks reject negative ages or
cell counts. Blank responses are stored as SQL `NULL`.

The design scales because new populations are added as rows rather than new
database columns. Indexes support the main filters for condition, treatment,
response, sample type, timepoint, and population. For substantially larger
datasets, the same schema can move to PostgreSQL and add batch loading,
partitioning, materialized summaries, or more indexes without changing the
logical model.

## Code structure

| File | Purpose |
| --- | --- |
| `load_data.py` | Create the SQLite schema and load `cell-count.csv` |
| `analyze_frequencies.py` | Produce the Part 2 frequency table |
| `analyze_response.py` | Run Part 3 tests and create boxplots |
| `analyze_subsets.py` | Run the Part 4 database queries |
| `dashboard.py` | Interactive Streamlit dashboard |
| `tests/` | Automated tests using temporary databases |
| `Makefile` | Codespaces setup, pipeline, dashboard, and test commands |

The scripts separate querying, calculation, output writing, and command-line
execution. Optional input and output paths allow tests to run without changing
the repository artifacts.

## Tests

```bash
make test
```

Tests cover data loading, schema integrity, validation, frequency calculations,
statistics, subset queries, output files, and dashboard rendering.

## Dashboard

Run:

```bash
make dashboard
```

The dashboard provides:

- trial overview charts;
- filters and per-sample relative frequencies;
- interactive miraclib response boxplots and statistics;
- baseline cohort counts and sample data.

Codespaces displays a forwarded dashboard URL for port 8501 after the server
starts. Open that link when prompted.
