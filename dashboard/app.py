# dashboard/app.py
# ─────────────────────────────────────────────────────────────
# Streamlit dashboard for genomic variant data (Y chromosome)
# Uses cleaned and aggregated Parquet files from data/processed/
# ─────────────────────────────────────────────────────────────

import streamlit as st
import pandas as pd
import plotly.express as px

# ==== PAGE SETUP ====
st.set_page_config(page_title="Genomic Variant Dashboard", layout="wide")
st.title("🧬 Genomic Variant Dashboard 🧬")

st.markdown("""
This dashboard uses **genomic variant data** to demonstrate:

- Population allele frequency maps
- Allele frequency distributions
- Variant density along a chromosome
- Interactive variant table
""")

# ==== SIDEBAR CONTROLS ====
st.sidebar.header("Controls")

# Population filter
pop_filter = st.sidebar.multiselect(
    "Select Populations",
    ["African", "American", "East Asian", "European", "South Asian"],
    default=["African", "European", "East Asian"]
)

# ==== LOAD DATA ====
@st.cache_data
def load_data():
    try:
        # Full cleaned variants table
        variants = pd.read_parquet("../data/processed/variants_cleaned.parquet")
        # Aggregated variant density per window
        density = pd.read_parquet("../data/processed/variants_aggregated.parquet")
        return variants, density
    except Exception as e:
        st.error(f"Error loading processed data: {e}")
        return None, None

variants, density = load_data()
if variants is None or density is None:
    st.stop()

# ==== MELT POPULATION AF FOR PLOTTING ====
pop_cols = ["af_afr", "af_amr", "af_eas", "af_eur", "af_sas"]
pop_labels = {
    "af_afr": "African",
    "af_amr": "American",
    "af_eas": "East Asian",
    "af_eur": "European",
    "af_sas": "South Asian"
}

df_melted = variants[["pos", "variant_type"] + pop_cols].melt(
    id_vars=["pos", "variant_type"],
    value_vars=pop_cols,
    var_name="population",
    value_name="allele_frequency"
)
df_melted["population"] = df_melted["population"].map(pop_labels)
df_melted = df_melted[df_melted["population"].isin(pop_filter)]

# ==== TABS ====
tab1, tab2, tab3, tab4 = st.tabs([
    "📊 Summary",
    "🌍 Global Map",
    "📈 Allele Frequency",
    "📍 Variant Density"
])

# ---- Summary ----
with tab1:
    st.header("Summary Metrics")
    col1, col2, col3 = st.columns(3)
    col1.metric("Total Variants", len(variants))
    col2.metric("SNP %", f"{(variants['variant_type']=='SNP').mean()*100:.1f}%")
    col3.metric("Mean AF", f"{df_melted['allele_frequency'].mean():.4f}")

# ---- Global Map ----
with tab2:
    st.header("🌍 Population Allele Frequency Map")
    coords = {"African": (0, 10), "American": (-60, 10), "East Asian": (105, 30),
              "European": (10, 50), "South Asian": (78, 22)}
    map_df = df_melted.groupby("population")["allele_frequency"].mean().reset_index()
    map_df["lon"] = map_df["population"].apply(lambda p: coords[p][0])
    map_df["lat"] = map_df["population"].apply(lambda p: coords[p][1])

    fig_map = px.scatter_geo(
        map_df, lon="lon", lat="lat",
        size="allele_frequency",
        color="population",
        hover_name="population",
        projection="natural earth",
        title="Population Mean Allele Frequencies"
    )
    st.plotly_chart(fig_map, use_container_width=True)

# ---- Allele Frequency ----
with tab3:
    st.header("Allele Frequency Distribution")
    fig_violin = px.violin(
        df_melted, x="population", y="allele_frequency",
        color="population", box=True, points="all",
        title="Allele Frequency by Population"
    )
    st.plotly_chart(fig_violin, use_container_width=True)

# ---- Variant Density ----
with tab4:
    st.header("Variant Density Across Chromosome")
    fig_density = px.bar(
        density,
        x=density["window"]/1_000_000,  # in Mb
        y="variant_count",
        labels={"x":"Position (Mb)", "variant_count":"Variant Count"},
        title="Variant Density per 1Mb Window",
        color="variant_count",
        color_continuous_scale="viridis"
    )
    st.plotly_chart(fig_density, use_container_width=True)

# ---- Variant Table ----
st.markdown("### 🧾 Variant Table (first 100 rows)")
st.dataframe(variants.head(100))
