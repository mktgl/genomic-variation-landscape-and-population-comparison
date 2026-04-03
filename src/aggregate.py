# src/aggregate.py
# ─────────────────────────────────────────────────────────────
# PURPOSE: Functions to summarise variants by genomic window and population
# Used by: notebooks/03_aggregation.ipynb
# ─────────────────────────────────────────────────────────────

import pandas as pd


def compute_window_density(df: pd.DataFrame, window_size: int = 100_000) -> pd.DataFrame:
    """
    Divides the chromosome into windows and counts variants per window.

    Parameters:
        df          : filtered variant DataFrame (must have 'pos', 'variant_type', 'af_global')
        window_size : size of each window in base pairs (default 100kb)

    Returns:
        DataFrame with one row per window: window, variant_count, mean_af, snp_count, indel_count
    """
    df = df.copy()
    df['window'] = (df['pos'] // window_size) * window_size

    result = df.groupby('window').agg(
        variant_count = ('pos', 'count'),
        mean_af       = ('af_global', 'mean'),
        snp_count     = ('variant_type', lambda x: (x == 'SNP').sum()),
        indel_count   = ('variant_type', lambda x: (x == 'INDEL').sum())
    ).reset_index()

    return result


def melt_population_af(df: pd.DataFrame) -> pd.DataFrame:
    """
    Reshapes the DataFrame from wide (one column per population) to long format
    (one row per variant-population pair). Useful for violin plots and comparisons.

    Returns:
        DataFrame with columns: pos, variant_type, population, allele_frequency
    """
    pop_cols = ['af_afr', 'af_eur', 'af_eas', 'af_amr', 'af_sas']
    pop_labels = {
        'af_afr': 'African',
        'af_eur': 'European',
        'af_eas': 'East Asian',
        'af_amr': 'American',
        'af_sas': 'South Asian'
    }

    result = df[['pos', 'variant_type'] + pop_cols].melt(
        id_vars=['pos', 'variant_type'],
        value_vars=pop_cols,
        var_name='population',
        value_name='allele_frequency'
    )
    result['population'] = result['population'].map(pop_labels)
    return result.dropna(subset=['allele_frequency'])
