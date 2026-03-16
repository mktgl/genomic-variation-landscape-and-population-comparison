# src/filter.py
# ─────────────────────────────────────────────────────────────
# PURPOSE: Functions to clean and filter genomic variant tables
# Used by: notebooks/02_cleaning_filtering.ipynb
# ─────────────────────────────────────────────────────────────

import pandas as pd


def filter_pass_only(df: pd.DataFrame) -> pd.DataFrame:
    """Keep only variants that passed gnomAD quality control (FILTER == PASS)."""
    return df[df['filter_status'] == 'PASS'].copy()


def drop_missing_af(df: pd.DataFrame) -> pd.DataFrame:
    """Remove rows where global allele frequency is missing."""
    return df.dropna(subset=['af_global']).copy()


def apply_maf_threshold(df: pd.DataFrame, min_maf: float = 0.001) -> pd.DataFrame:
    """
    Keep only variants with allele frequency >= min_maf.

    Parameters:
        df      : input DataFrame (must have 'af_global' column)
        min_maf : minimum allele frequency to keep (default 0.001 = 0.1%)

    Returns:
        Filtered DataFrame
    """
    return df[df['af_global'] >= min_maf].copy()


def full_filter_pipeline(df: pd.DataFrame, min_maf: float = 0.001) -> pd.DataFrame:
    """
    Runs all three filters in order:
      1. PASS only
      2. Drop missing AF
      3. MAF threshold

    Returns cleaned DataFrame.
    """
    before = len(df)
    df = filter_pass_only(df)
    df = drop_missing_af(df)
    df = apply_maf_threshold(df, min_maf)
    print(f'Filtering: {before} → {len(df)} variants (removed {before - len(df)})')
    return df
