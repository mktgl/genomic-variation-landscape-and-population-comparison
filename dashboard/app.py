import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px

# ==== PAGE SETUP ====
st.set_page_config(page_title="Genomic Variant Dashboard", layout="wide")
st.title("🧬 Genomic Variant Dashboard (Synthetic Data)")

st.markdown("""
This dashboard uses **synthetic genomic variant data** to demonstrate:
- Population allele frequency maps
- Allele frequency distributions
- Variant density along a chromosome
- Interactive variant table
""")

# ==== SIDEBAR CONTROLS ====
st.sidebar.header("Controls")
num_variants = st.sidebar.slider("Number of synthetic variants", 100, 5000, 1000)
pop_filter = st.sidebar.multiselect(
    "Select Populations",
    ["AFR", "AMR", "EAS", "EUR", "SAS"],
    default=["AFR", "EUR", "EAS"]
)

# ==== GENERATE SYNTHETIC DATA ====
np.random.seed(42)

# Chromosome positions
positions = np.random.randint(1, 100_000_000, num_variants)
variant_types = np.random.choice(["SNP", "INDEL"], num_variants, p=[0.8, 0.2])

# Allele frequencies for populations
allele_freqs = {pop: np.random.beta(a=1, b=10, size=num_variants) for pop in ["AFR", "AMR", "EAS", "EUR", "SAS"]}

df = pd.DataFrame({
    "position": positions,
    "variant_type": variant_types,
    **allele_freqs
})

# Ensure numeric types
df["position"] = df["position"].astype(int)

# Melt for plotting
df_melted = df.melt(id_vars=["position", "variant_type"], 
                    value_vars=["AFR","AMR","EAS","EUR","SAS"], 
                    var_name="population", value_name="allele_frequency")
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
    col1.metric("Total Variants", len(df))
    col2.metric("SNP %", f"{(df['variant_type']=='SNP').mean()*100:.1f}%")
    col3.metric("Mean AF", f"{df_melted['allele_frequency'].mean():.4f}")

# ---- Global Map ----
with tab2:
    st.header("🌍 Population Allele Frequency Map")
    coords = {"AFR": (0, 10), "AMR": (-60, 10), "EAS": (105, 30), "EUR": (10, 50), "SAS": (78, 22)}
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
    # Bin positions in 1Mb windows
    df["window"] = (df["position"] // 1_000_000) * 1_000_000
    density = df.groupby("window").size().reset_index(name="variant_count")
    
    fig_density = px.bar(
        density,
        x=density["window"]/1_000_000,  # divide numeric column
        y="variant_count",
        labels={"x":"Position (Mb)", "variant_count":"Variant Count"},
        title="Variant Density per 1Mb Window",
        color="variant_count",
        color_continuous_scale="viridis"
    )
    st.plotly_chart(fig_density, use_container_width=True)

# ---- Variant Table ----
st.markdown("### 🧾 Variant Table (first 100 rows)")
st.dataframe(df.head(100))
