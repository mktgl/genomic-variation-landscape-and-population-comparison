the interactive web dashboard for the project

# ─────────────────────────────────────────────────────────────
# PURPOSE: Interactive web dashboard built with Streamlit
# Run with: streamlit run dashboard/app.py
#
# This is the FINAL product of the project — an interactive
# webpage where anyone can explore the chrY variant data
# without writing any code.
# ─────────────────────────────────────────────────────────────

import streamlit as st
import pandas as pd
import plotly.express as px
import os

st.set_page_config(page_title='Genomic Variant Dashboard', layout='wide')

st.title('🧬 Genomic Variant Landscape — chrY (gnomAD v4.1)')
st.markdown('Explore variant density and population allele frequencies along the human Y chromosome.')

# ── LOAD DATA ──────────────────────────────────────────────────
@st.cache_data
def load_data():
    base = os.path.join(os.path.dirname(__file__), '..', 'data', 'processed')
    window_df = pd.read_parquet(os.path.join(base, 'window_density.parquet'))
    pop_df    = pd.read_parquet(os.path.join(base, 'population_af.parquet'))
    df        = pd.read_parquet(os.path.join(base, 'chrY_filtered.parquet'))
    return window_df, pop_df, df

try:
    window_df, pop_df, df = load_data()
    data_loaded = True
except FileNotFoundError:
    st.warning('⚠️  Processed data not found. Please run notebooks 01–03 first.')
    data_loaded = False

if data_loaded:
    # ── SIDEBAR FILTERS ────────────────────────────────────────
    st.sidebar.header('Filters')
    maf_min = st.sidebar.slider('Minimum Allele Frequency (MAF)', 0.0, 1.0, 0.001, 0.001)
    vtype   = st.sidebar.multiselect('Variant Type', ['SNP', 'INDEL'], default=['SNP', 'INDEL'])

    filtered = df[(df['af_global'] >= maf_min) & (df['variant_type'].isin(vtype))]
    st.sidebar.metric('Variants shown', f'{len(filtered):,}')

    # ── CHART 1: Density ───────────────────────────────────────
    st.subheader('📍 Variant Density Along chrY')
    fig1 = px.bar(window_df, x=window_df['window']/1e6, y='variant_count',
                  labels={'x': 'Genomic Position (Mb)', 'variant_count': 'Variants per 100kb'},
                  color='mean_af', color_continuous_scale='Blues')
    st.plotly_chart(fig1, use_container_width=True)

    # ── CHART 2: Population AF ─────────────────────────────────
    st.subheader('🌍 Allele Frequency by Population')
    fig2 = px.violin(pop_df, x='population', y='allele_frequency', box=True,
                     color='population',
                     labels={'allele_frequency': 'Allele Frequency', 'population': 'Population'})
    st.plotly_chart(fig2, use_container_width=True)

    # ── CHART 3: AFR vs EUR ────────────────────────────────────
    st.subheader('🔬 Population Comparison: AFR vs EUR')
    scatter_df = filtered[['af_afr', 'af_eur', 'variant_type']].dropna()
    fig3 = px.scatter(scatter_df, x='af_afr', y='af_eur', color='variant_type',
                      opacity=0.4,
                      labels={'af_afr': 'AF (African)', 'af_eur': 'AF (European)'})
    st.plotly_chart(fig3, use_container_width=True)

    # ── RAW TABLE ──────────────────────────────────────────────
    with st.expander('View raw data table'):
        st.dataframe(filtered.head(200))
