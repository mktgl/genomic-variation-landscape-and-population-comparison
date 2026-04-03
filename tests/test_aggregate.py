# tests/test_aggregate.py
import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
import pandas as pd
from src.aggregate import compute_window_density, melt_population_af

@pytest.fixture
def sample_variant_df():
    """Sample DataFrame for testing aggregation functions"""
    data = {
        'pos': [100_000, 150_000, 250_000, 300_000, 450_000],
        'variant_type': ['SNP', 'INDEL', 'SNP', 'SNP', 'INDEL'],
        'af_global': [0.05, 0.12, 0.003, 0.08, 0.25],
        'af_afr': [0.1, 0.15, 0.01, 0.09, 0.3],
        'af_eur': [0.02, 0.1, 0.0, 0.07, 0.2],
        'af_eas': [0.03, 0.08, 0.005, 0.06, 0.22],
        'af_amr': [0.04, 0.14, None, 0.085, 0.28],
        'af_sas': [0.06, 0.11, 0.002, 0.075, 0.24],
    }
    return pd.DataFrame(data)


def test_compute_window_density_default_window(sample_variant_df):
    result = compute_window_density(sample_variant_df, window_size=100_000)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 3  # windows: 100k, 200k, 400k
    assert list(result.columns) == ['window', 'variant_count', 'mean_af', 'snp_count', 'indel_count']
    
    # Check window calculations
    assert result.iloc[0]['window'] == 100_000
    assert result.iloc[0]['variant_count'] == 2
    assert result.iloc[0]['snp_count'] == 1
    assert result.iloc[0]['indel_count'] == 1
    
    assert result.iloc[2]['window'] == 400_000
    assert result.iloc[2]['variant_count'] == 1


def test_compute_window_density_custom_window(sample_variant_df):
    result = compute_window_density(sample_variant_df, window_size=200_000)
    
    assert len(result) == 2  # windows: 0 and 400k
    assert result.iloc[0]['variant_count'] == 3
    assert result.iloc[1]['variant_count'] == 1


def test_melt_population_af(sample_variant_df):
    result = melt_population_af(sample_variant_df)
    
    assert isinstance(result, pd.DataFrame)
    assert len(result) == 25  # 5 variants × 5 populations
    
    assert list(result.columns) == ['pos', 'variant_type', 'population', 'allele_frequency']
    
    # Check population label mapping
    assert set(result['population'].unique()) == {
        'African', 'European', 'East Asian', 'American', 'South Asian'
    }
    
    # Check no NaN allele frequencies
    assert result['allele_frequency'].isna().sum() == 0


def test_melt_population_af_drops_na(sample_variant_df):
    # Introduce a NaN in one population
    df_with_na = sample_variant_df.copy()
    df_with_na.loc[0, 'af_amr'] = None
    
    result = melt_population_af(df_with_na)
    
    # Should drop the NaN entry
    assert len(result) == 24  # 25 - 1
