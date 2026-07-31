"""Interactive Streamlit dashboard for the Loblaw Bio analysis."""

from __future__ import annotations

import sqlite3
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

from analyze_response import (
    SUBJECT_FREQUENCY_QUERY,
    calculate_statistics,
)
from analyze_subsets import (
    BASELINE_SAMPLES_QUERY,
    B_CELL_AVERAGE_QUERY,
    GENDER_COUNTS_QUERY,
    PROJECT_COUNTS_QUERY,
    RESPONSE_COUNTS_QUERY,
)
from load_data import POPULATIONS


ROOT = Path(__file__).resolve().parent
DB_PATH = ROOT / "cell-count.db"

st.set_page_config(
    page_title="Loblaw Bio | Immune Cell Analysis",
    page_icon="🧬",
    layout="wide",
)


@st.cache_data
def query_dataframe(query: str) -> pd.DataFrame:
    """Run a read-only query against the generated SQLite database."""
    with sqlite3.connect(DB_PATH) as connection:
        return pd.read_sql_query(query, connection)


@st.cache_data
def load_frequency_data() -> pd.DataFrame:
    return query_dataframe(
        """
        WITH sample_totals AS (
            SELECT sample_id, SUM(cell_count) AS total_count
            FROM cell_counts
            GROUP BY sample_id
        )
        SELECT
            projects.project_name AS project,
            subjects.subject_code AS subject,
            subjects.condition,
            subjects.gender,
            subjects.treatment,
            COALESCE(subjects.response, 'not recorded') AS response,
            samples.sample_code AS sample,
            samples.sample_type,
            samples.time_from_treatment_start AS timepoint,
            sample_totals.total_count,
            cell_populations.population_name AS population,
            cell_counts.cell_count AS count,
            CASE
                WHEN sample_totals.total_count = 0 THEN 0.0
                ELSE 100.0 * cell_counts.cell_count / sample_totals.total_count
            END AS percentage
        FROM samples
        JOIN subjects USING (subject_id)
        JOIN projects USING (project_id)
        JOIN sample_totals USING (sample_id)
        JOIN cell_counts USING (sample_id)
        JOIN cell_populations USING (population_id)
        ORDER BY samples.sample_code, cell_populations.population_id
        """
    )


@st.cache_data
def load_response_data() -> tuple[pd.DataFrame, pd.DataFrame]:
    subject_frequencies = query_dataframe(SUBJECT_FREQUENCY_QUERY)
    groups = {
        population: {"yes": [], "no": []} for population in POPULATIONS
    }
    for row in subject_frequencies.itertuples(index=False):
        groups[row.population][row.response].append(float(row.mean_percentage))
    statistics = pd.DataFrame(calculate_statistics(groups))
    return subject_frequencies, statistics


def multiselect_filter(
    frame: pd.DataFrame, column: str, label: str
) -> pd.DataFrame:
    options = sorted(frame[column].dropna().unique().tolist())
    selected = st.multiselect(label, options, default=options)
    return frame[frame[column].isin(selected)]


if not DB_PATH.is_file():
    st.error("Database not found. Run `make pipeline` before starting the dashboard.")
    st.stop()

st.title("Loblaw Bio immune cell analysis")
st.caption(
    "Explore sample composition, miraclib response patterns, and baseline "
    "melanoma cohorts."
)

frequencies = load_frequency_data()
response_data, response_statistics = load_response_data()
baseline_samples = query_dataframe(BASELINE_SAMPLES_QUERY)

overview_tab, frequency_tab, response_tab, baseline_tab = st.tabs(
    (
        "Trial overview",
        "Cell frequencies",
        "Miraclib response",
        "Baseline subset",
    )
)

with overview_tab:
    metric_columns = st.columns(4)
    metric_columns[0].metric("Projects", frequencies["project"].nunique())
    metric_columns[1].metric("Subjects", frequencies["subject"].nunique())
    metric_columns[2].metric("Samples", frequencies["sample"].nunique())
    metric_columns[3].metric("Cell populations", frequencies["population"].nunique())

    sample_overview = frequencies.drop_duplicates("sample")
    left, right = st.columns(2)
    with left:
        project_counts = (
            sample_overview.groupby("project", as_index=False)
            .size()
            .rename(columns={"size": "samples"})
        )
        figure = px.bar(
            project_counts,
            x="project",
            y="samples",
            text_auto=True,
            title="Samples by project",
            labels={"project": "Project", "samples": "Samples"},
        )
        figure.update_layout(showlegend=False)
        st.plotly_chart(figure, width="stretch")
    with right:
        condition_counts = (
            sample_overview.groupby("condition", as_index=False)
            .size()
            .rename(columns={"size": "samples"})
        )
        figure = px.bar(
            condition_counts,
            x="condition",
            y="samples",
            text_auto=True,
            title="Samples by condition",
            labels={"condition": "Condition", "samples": "Samples"},
        )
        figure.update_layout(showlegend=False)
        st.plotly_chart(figure, width="stretch")

with frequency_tab:
    st.subheader("Relative frequency by sample")
    st.caption(
        "Filter the trial metadata, then select a sample to inspect the five "
        "cell populations."
    )
    filter_columns = st.columns(4)
    filtered = frequencies
    with filter_columns[0]:
        filtered = multiselect_filter(filtered, "project", "Project")
    with filter_columns[1]:
        filtered = multiselect_filter(filtered, "condition", "Condition")
    with filter_columns[2]:
        filtered = multiselect_filter(filtered, "treatment", "Treatment")
    with filter_columns[3]:
        filtered = multiselect_filter(filtered, "sample_type", "Sample type")

    available_samples = sorted(filtered["sample"].unique().tolist())
    if not available_samples:
        st.warning("No samples match the selected filters.")
    else:
        selected_sample = st.selectbox("Sample", available_samples)
        selected = filtered[filtered["sample"] == selected_sample].copy()
        chart, table = st.columns((3, 2))
        with chart:
            figure = px.bar(
                selected,
                x="population",
                y="percentage",
                text="percentage",
                labels={
                    "population": "Cell population",
                    "percentage": "Relative frequency (%)",
                },
            )
            figure.update_traces(texttemplate="%{text:.2f}%", textposition="outside")
            figure.update_layout(showlegend=False)
            st.plotly_chart(figure, width="stretch")
        with table:
            st.dataframe(
                selected[
                    ["population", "count", "total_count", "percentage"]
                ].style.format({"percentage": "{:.2f}%"}),
                hide_index=True,
                width="stretch",
            )

with response_tab:
    st.subheader("Melanoma PBMC frequencies by miraclib response")
    st.caption(
        "Each point is a subject's mean frequency across eligible PBMC "
        "timepoints. Tests are two-sided Mann–Whitney U with "
        "Benjamini–Hochberg FDR correction."
    )
    plot_data = response_data.copy()
    plot_data["Response"] = plot_data["response"].map(
        {"yes": "Responders", "no": "Non-responders"}
    )
    plot_data["Population"] = (
        plot_data["population"].str.replace("_", " ").str.title()
    )
    figure = px.box(
        plot_data,
        x="Population",
        y="mean_percentage",
        color="Response",
        points="outliers",
        labels={"mean_percentage": "Subject mean relative frequency (%)"},
        category_orders={"Response": ["Responders", "Non-responders"]},
    )
    st.plotly_chart(figure, width="stretch")

    display_statistics = response_statistics[
        [
            "population",
            "responder_subjects",
            "non_responder_subjects",
            "responder_median_percentage",
            "non_responder_median_percentage",
            "p_value",
            "adjusted_p_value",
            "significant",
            "rank_biserial_correlation",
        ]
    ].copy()
    display_statistics.columns = [
        "Population",
        "Responders",
        "Non-responders",
        "Responder median (%)",
        "Non-responder median (%)",
        "p-value",
        "FDR q-value",
        "Significant",
        "Rank-biserial effect",
    ]
    st.dataframe(
        display_statistics.style.format(
            {
                "Responder median (%)": "{:.2f}",
                "Non-responder median (%)": "{:.2f}",
                "p-value": "{:.4g}",
                "FDR q-value": "{:.4g}",
                "Rank-biserial effect": "{:.3f}",
            }
        ),
        hide_index=True,
        width="stretch",
    )
    st.info(
        "No population meets the FDR-adjusted significance threshold of 0.05. "
        "CD4 T cells have the smallest raw p-value, but q = 0.0621."
    )

with baseline_tab:
    st.subheader("Baseline melanoma, miraclib, PBMC subset")
    project_counts = dict(
        query_dataframe(PROJECT_COUNTS_QUERY).itertuples(index=False, name=None)
    )
    response_counts = dict(
        query_dataframe(RESPONSE_COUNTS_QUERY).itertuples(index=False, name=None)
    )
    gender_counts = dict(
        query_dataframe(GENDER_COUNTS_QUERY).itertuples(index=False, name=None)
    )
    b_cell_average = float(query_dataframe(B_CELL_AVERAGE_QUERY).iloc[0, 0])

    metrics = st.columns(4)
    metrics[0].metric("Baseline samples", len(baseline_samples))
    metrics[1].metric("Responders", response_counts.get("yes", 0))
    metrics[2].metric("Non-responders", response_counts.get("no", 0))
    metrics[3].metric("Male responder B-cell average", f"{b_cell_average:.2f}")

    left, right = st.columns(2)
    with left:
        project_frame = pd.DataFrame(
            project_counts.items(), columns=["Project", "Samples"]
        )
        figure = px.bar(
            project_frame,
            x="Project",
            y="Samples",
            text_auto=True,
            title="Baseline samples by project",
        )
        figure.update_layout(showlegend=False)
        st.plotly_chart(figure, width="stretch")
    with right:
        subject_frame = pd.DataFrame(
            [
                ("Responders", response_counts.get("yes", 0)),
                ("Non-responders", response_counts.get("no", 0)),
                ("Female", gender_counts.get("F", 0)),
                ("Male", gender_counts.get("M", 0)),
            ],
            columns=["Group", "Subjects"],
        )
        figure = px.bar(
            subject_frame,
            x="Group",
            y="Subjects",
            text_auto=True,
            title="Baseline distinct-subject counts",
        )
        figure.update_layout(showlegend=False)
        st.plotly_chart(figure, width="stretch")

    st.dataframe(baseline_samples, hide_index=True, width="stretch")
