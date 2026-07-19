"""Rejuvenation Explorer — interactive companion app.

Lets wet-lab and computational collaborators explore every module's outputs
without touching code. Run:

    streamlit run app/streamlit_app.py -- --results results/
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import plotly.express as px
import streamlit as st

RESULTS = Path("results")
if "--results" in sys.argv:
    RESULTS = Path(sys.argv[sys.argv.index("--results") + 1])

st.set_page_config(page_title="Rejuvenation Explorer", layout="wide")
st.title("Cell Rejuvenation Atlas — Explorer")
st.caption(
    "Interactive views over reproducible reanalyses of peer-reviewed aging & "
    "reprogramming data (see DATASETS.md in the repo)."
)

module = st.sidebar.radio(
    "Module",
    [
        "1 · Epigenetic clock (Gill 2022)",
        "2 · Reprogramming fates (Schiebinger 2019)",
        "3 · Aging vs. reprogramming (Browder 2022 × TMS)",
        "4 · Spatial brain aging (Allen 2023)",
        "5 · SASP proteomics (Basisty 2020)",
    ],
)


@st.cache_data
def load_csv(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path) if path.exists() else None


def missing(name: str) -> None:
    st.info(
        f"No results found for {name} yet. Run the pipeline first:\n\n"
        f"`nextflow run main.nf -profile docker --module <name>`",
        icon="ℹ️",
    )


if module.startswith("1"):
    st.header("Epigenetic age across transient reprogramming")
    df = load_csv(RESULTS / "module1" / "dnam_age_timecourse.csv")
    stats_path = RESULTS / "module1" / "delta_age_stats.json"
    if df is None:
        missing("rejuvenation_clock")
    else:
        fig = px.line(df, x="day", y="dnam_age", markers=True,
                      labels={"day": "MPTR day", "dnam_age": "DNAm age (years)"},
                      title="GSE165179 — DNAm age falls below baseline after MPTR")
        st.plotly_chart(fig, use_container_width=True)
        if stats_path.exists():
            stats = json.loads(stats_path.read_text())
            cols = st.columns(3)
            cols[0].metric("Δ age (peak vs day 0)", f"{stats.get('delta_years', 0):.1f} y")
            cols[1].metric("95% CI", f"[{stats.get('ci95', [0, 0])[0]:.1f}, {stats.get('ci95', [0, 0])[1]:.1f}]")
            cols[2].metric("Cohen's dz", f"{stats.get('cohens_dz', 0):.2f}")

elif module.startswith("2"):
    st.header("Fate probabilities from optimal transport")
    df = load_csv(RESULTS / "module2" / "fate_probabilities.csv")
    if df is None:
        missing("trajectory")
    else:
        fig = px.histogram(df, x="ipsc_fate_prob", nbins=50,
                           labels={"ipsc_fate_prob": "P(reach iPSC fate)"},
                           title="Day-0 cells: probability of successful reprogramming")
        st.plotly_chart(fig, use_container_width=True)
        st.image(str(RESULTS / "module2" / "fate_analysis.png")
                 if (RESULTS / "module2" / "fate_analysis.png").exists() else None)

elif module.startswith("3"):
    st.header("Does partial reprogramming move aged tissue toward a young state?")
    df = load_csv(RESULTS / "module3" / "tissue_scores.csv")
    if df is None:
        missing("integration")
    else:
        st.caption("Browder 2022 bulk RNA-seq (6 tissues, 7-month cyclic OSK) projected "
                   "onto the Tabula Muris Senis aging axis. Negative Δ = 4F shifted young.")
        st.dataframe(df, use_container_width=True)
        png = RESULTS / "module3" / "aging_score.png"
        if png.exists():
            st.image(str(png))

elif module.startswith("4"):
    st.header("Spatial neighborhood changes with age (MERFISH)")
    files = sorted((RESULTS / "module4").glob("nhood_enrichment_*.csv")) if (RESULTS / "module4").exists() else []
    if not files:
        missing("spatial")
    else:
        age = st.selectbox("Age group", [f.stem.replace("nhood_enrichment_", "") for f in files])
        z = pd.read_csv(RESULTS / "module4" / f"nhood_enrichment_{age}.csv", index_col=0)
        fig = px.imshow(z, color_continuous_scale="RdBu_r", zmin=-z.abs().max().max(),
                        zmax=z.abs().max().max(), title=f"Neighborhood enrichment z-scores — {age}")
        st.plotly_chart(fig, use_container_width=True)

else:
    st.header("Core SASP from proteomics")
    df = load_csv(RESULTS / "module5" / "core_sasp.csv")
    if df is None:
        missing("sasp")
    else:
        st.metric("Core SASP proteins", len(df))
        st.dataframe(df, use_container_width=True)
        png = RESULTS / "module5" / "sasp_per_inducer.png"
        if png.exists():
            st.image(str(png))

st.sidebar.markdown("---")
st.sidebar.caption("Built with Streamlit · figures regenerated from code, never pasted")
