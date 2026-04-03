import streamlit as st
import pandas as pd
import plotly.express as px
import sys
import os

# ====================== PATH SETUP ======================
sys.path.insert(0, os.path.abspath('..'))   # So we can import from src/

from src.ingest import parse_vcf_to_dataframe
from src.filter import full_filter_pipeline
from src.aggregate import compute_window_density, melt_population_af

# ====================== PAGE CONFIG ======================
st.set_page_config(page_title="Genomic Variant Dashboard", layout="wide")
st.title("🧬 Genomic Variant Landscape & Population Comparison")
st.markdown("**Live analysis** using your `src/` modules | chrY gnomAD v4.1")

# ====================== SIDEBAR ======================
st.sidebar.header("🎛️ Controls")

uploaded_file = st.sidebar.file_uploader(
    "Upload VCF file (.vcf or .vcf.gz)", 
    type=["vcf", "vcf.gz"]
)

max_rows = st.sidebar.slider("Max variants to load (for speed)", 
                             min_value=500, max_value=20000, value=5000, step=500)

min_maf = st.sidebar.slider("Minimum Global MAF", 
                            min_value=0.0, max_value=0.1, value=0.001, step=0.0001)

variant_types = st.sidebar.multiselect(
    "Variant Type", 
    options=["SNP", "INDEL"], 
    default=["SNP", "INDEL"]
)

# ====================== MAIN DASHBOARD ======================
if uploaded_file is None:
    st.info("👆 Upload a VCF file to start analysis.")
    st.stop()

# Save uploaded file temporarily
temp_path = "temp_uploaded.vcf"
with open(temp_path, "wb") as f:
    f.write(uploaded_file.getbuffer())

# ------------------- Ingest -------------------
with st.spinner("Parsing VCF file..."):
    df_raw = parse_vcf_to_dataframe(temp_path, max_rows=max_rows)

st.success(f"✅ Loaded **{len(df_raw):,}** variants from VCF")

# ------------------- Filter -------------------
with st.spinner("Applying filters..."):
    df_filtered = full_filter_pipeline(df_raw, min_maf=min_maf)

# Additional variant type filter
df_filtered = df_filtered[df_filtered['variant_type'].isin(variant_types)]

st.metric("Variants after filtering", f"{len(df_filtered):,}")

if df_filtered.empty:
    st.error("No variants remaining after filtering. Try lowering MAF or removing variant type filter.")
    st.stop()

# ------------------- Tabs -------------------
tab1, tab2, tab3, tab4 = st.tabs([
    "📊 Overview", 
    "🌍 Population Allele Frequency", 
    "🪟 Window Density", 
    "📋 Raw Data"
])

with tab1:
    st.subheader("Summary Statistics")
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Total Raw", len(df_raw))
    col2.metric("After Filter", len(df_filtered))
    col3.metric("SNP %", f"{(df_filtered['variant_type']=='SNP').mean()*100:.1f}%")
    col4.metric("Mean AF", f"{df_filtered['af_global'].mean():.4f}")

    st.dataframe(df_filtered.head(10), use_container_width=True)

with tab2:
    st.subheader("Allele Frequency Distribution by Population")
    melted = melt_population_af(df_filtered)
    
    fig = px.violin(melted, x='population', y='allele_frequency', 
                    color='population', box=True, points="outliers",
                    title="Population-wise Allele Frequency")
    st.plotly_chart(fig, use_container_width=True)

    # AFR vs EUR scatter
    st.subheader("African vs European AF Comparison")
    scatter = df_filtered[['af_afr', 'af_eur', 'variant_type']].dropna()
    fig2 = px.scatter(scatter, x='af_afr', y='af_eur', color='variant_type',
                      opacity=0.5, title="AF African vs European")
    st.plotly_chart(fig2, use_container_width=True)

with tab3:
    st.subheader("Variant Density by Genomic Window (100kb)")
    agg_df = compute_window_density(df_filtered, window_size=100_000)
    
    fig3 = px.bar(agg_df, x=agg_df['window']/1_000_000, y='variant_count',
                  color='mean_af', color_continuous_scale='Blues',
                  labels={'x': 'Position on chrY (Mb)', 'variant_count': 'Number of Variants'})
    st.plotly_chart(fig3, use_container_width=True)
    
    st.dataframe(agg_df.head(15), use_container_width=True)

with tab4:
    st.subheader("Filtered Variants Table")
    st.dataframe(df_filtered, use_container_width=True)

# Download button
csv = df_filtered.to_csv(index=False).encode()
st.download_button("📥 Download Filtered Data as CSV", 
                   csv, "filtered_variants.csv", "text/csv")

# Cleanup temp file
if os.path.exists(temp_path):
    os.remove(temp_path)
