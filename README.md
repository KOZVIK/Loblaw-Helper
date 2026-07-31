# Loblaw Bio Immune Cell Analysis

This project loads clinical-trial immune cell counts into a normalized SQLite
database and provides reproducible analyses of cell-population frequencies,
miraclib treatment response, and baseline patient subsets.

## Requirements

- Python 3.10 or newer
- `pip`

The database loader and Parts 2 and 4 use only the Python standard library.
Part 3 requires SciPy and Matplotlib for statistical testing and plotting.

## Setup

In GitHub Codespaces or another environment with GNU Make:

```bash
make setup
```

In Windows PowerShell without GNU Make:

```powershell
python -m pip install --requirement requirements.txt
```

This runs `python -m pip install --requirement requirements.txt`. The
`--requirement` option tells pip to read the dependency list from the file; it
does not refer to the R programming language.

## Reproducing the analysis

Run the complete pipeline:

```bash
make pipeline
```

This executes `load_data.py`, `analyze_frequencies.py`, `analyze_response.py`,
and `analyze_subsets.py` sequentially. The scripts do not require command-line
arguments. `load_data.py` reads `cell-count.csv` relative to its own location,
so it can also be invoked from a different working directory.

In Windows PowerShell without GNU Make, run the same pipeline directly:

```powershell
python load_data.py
python analyze_frequencies.py
python analyze_response.py
python analyze_subsets.py
```

### Part 1: Load the database

```bash
python load_data.py
```

This recreates `cell-count.db` in the repository root and loads:

- 3 projects
- 3,500 subjects
- 10,500 samples
- 5 cell populations
- 52,500 cell-count records

Running the loader again safely rebuilds the generated tables without
duplicating records. Empty response values are represented as SQL `NULL`.

### Part 2: Cell-population frequencies

```bash
python analyze_frequencies.py
```

For every sample, the script sums the five population counts and calculates:

```text
percentage = 100 * population count / total sample count
```

It prints a preview and writes
[`outputs/cell_frequencies.csv`](outputs/cell_frequencies.csv). The output has
one row per sample and population and contains the required columns:

| Column | Description |
| --- | --- |
| `sample` | Source sample identifier |
| `total_count` | Sum of the five population counts |
| `population` | Immune cell population |
| `count` | Population cell count |
| `percentage` | Population percentage of the sample total |

### Part 3: Miraclib response analysis

```bash
python analyze_response.py
```

The analysis includes only melanoma subjects treated with miraclib whose
samples are PBMCs and whose response is `yes` or `no`.

Because each subject can have repeated longitudinal samples, the program first
calculates each subject's mean relative frequency across their eligible
timepoints. It then compares responders and non-responders using a two-sided
Mann–Whitney U test for each population. Benjamini–Hochberg false-discovery-rate
correction accounts for testing five populations. Statistical significance is
defined as an adjusted p-value below 0.05.

The analyzed groups contain 331 responders and 325 non-responders. No cell
population is significant after correction. CD4 T cells have the smallest raw
p-value (`p = 0.0124`), but the adjusted value is `q = 0.0621`.

Generated files:

- [`outputs/response_boxplots.png`](outputs/response_boxplots.png)
- [`outputs/response_statistics.csv`](outputs/response_statistics.csv)
- [`outputs/response_statistics.txt`](outputs/response_statistics.txt)

The statistics table reports group sizes, medians, median differences,
Mann–Whitney U statistics, raw and adjusted p-values, significance decisions,
and rank-biserial effect sizes.

### Part 4: Baseline subset analysis

```bash
python analyze_subsets.py
```

The baseline subset contains melanoma PBMC samples collected at time 0 from
subjects treated with miraclib:

| Measure | Result |
| --- | ---: |
| Total samples | 656 |
| `prj1` samples | 384 |
| `prj3` samples | 272 |
| Responder subjects | 331 |
| Non-responder subjects | 325 |
| Female subjects | 312 |
| Male subjects | 344 |

For responding melanoma males at time 0, across all sample and treatment types,
the average B-cell count is **10206.15**.

Generated files:

- [`outputs/baseline_melanoma_miraclib_pbmc.csv`](outputs/baseline_melanoma_miraclib_pbmc.csv)
- [`outputs/baseline_subset_counts.csv`](outputs/baseline_subset_counts.csv)
- [`outputs/baseline_male_responder_b_cell_average.txt`](outputs/baseline_male_responder_b_cell_average.txt)

## Database schema

The SQLite database uses five tables:

```text
projects 1 ──< subjects 1 ──< samples 1 ──< cell_counts >── 1 cell_populations
```

- `projects` stores each project once.
- `subjects` stores subject-level metadata: source subject code, condition, age,
  gender, treatment, and response. A subject belongs to one project.
- `samples` stores the sample identifier, sample type, and time from treatment
  start. Each sample belongs to one subject.
- `cell_populations` is a lookup table containing population names.
- `cell_counts` is the junction table between samples and populations and
  stores the observed count.

The schema preserves the source field name `condition`. Surrogate integer
primary keys provide efficient joins, while uniqueness constraints retain the
identity of source projects, subjects, samples, and populations. Foreign keys
enforce relationships, and non-negative checks protect ages and cell counts.

### Why this design scales

Project, subject, and sample metadata are not repeated for every population
measurement. A new immune population can be inserted into `cell_populations`
without altering the database schema, and a new measurement becomes one
`cell_counts` row. This is more maintainable than adding a database column for
every future population.

Indexes support the principal filtering patterns:

- condition, treatment, and response at the subject level;
- subject, timepoint, and sample type at the sample level;
- population-level cell-count queries.

For hundreds of projects and thousands or millions of samples, the same logical
model can be migrated to a server database such as PostgreSQL. Batch insertion,
table partitioning by project, materialized summary views, and additional
indexes can then be introduced based on observed query patterns without
redesigning the data model.

## Code structure

| File | Responsibility |
| --- | --- |
| `load_data.py` | Validate the CSV, initialize the schema, and load the database |
| `analyze_frequencies.py` | Generate the Part 2 sample-frequency table |
| `analyze_response.py` | Run Part 3 statistical tests and create boxplots |
| `analyze_subsets.py` | Run the Part 4 baseline subset queries |
| `dashboard.py` | Serve the interactive Streamlit dashboard |
| `Makefile` | Install dependencies and orchestrate the complete pipeline |
| `tests/` | Unit and integration tests using isolated temporary databases |
| `outputs/` | Reproducible tables, report text, and plots |

Each analysis module separates database queries, calculation, output writing,
and command-line orchestration. Functions accept optional database and output
paths so tests can run in temporary directories without modifying repository
artifacts.

## Tests

Run the complete test suite with:

```bash
make test
```

Windows PowerShell equivalent:

```powershell
python -m unittest discover -s tests -v
```

The suite covers schema integrity, normalized loading, null handling,
idempotency, input validation, frequency calculations, multiple-testing
correction, response-group validation, Part 4 filtering, grouped counts, and
the requested B-cell average.

## Dashboard

After running `make setup` and `make pipeline`, start the interactive dashboard:

```bash
make dashboard
```

Windows PowerShell equivalent:

```powershell
python -m streamlit run dashboard.py
```

Streamlit prints the local and forwarded URLs in the terminal. In GitHub
Codespaces, open the forwarded port when prompted. The dashboard includes:

- a trial overview by project and condition;
- interactive metadata filters and per-sample population frequencies;
- responder/non-responder boxplots and the complete statistical results;
- Part 4 baseline cohort metrics, charts, and sample-level data.

Dashboard link: run the command above to create the environment-specific local
or GitHub Codespaces URL.
