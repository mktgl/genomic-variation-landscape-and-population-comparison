# tests/test_aggregate.py

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
import pandas as pd
from src.aggregate import compute_window_density, melt_population_af


@pytest.fixture
def sample_variant_df():
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


# ─────────────────────────────────────────
# TEST: compute_window_density
# ─────────────────────────────────────────

def test_compute_window_density_default_window(sample_variant_df):
    result = compute_window_density(sample_variant_df, window_size=100_000)

    assert isinstance(result, pd.DataFrame)

    # Expected windows: 100k, 200k, 300k, 400k
    assert len(result) == 4
    assert list(result.columns) == ['window', 'variant_count', 'mean_af', 'snp_count', 'indel_count']

    # Check first window (100k)
    row = result[result['window'] == 100_000].iloc[0]
    assert row['variant_count'] == 2
    assert row['snp_count'] == 1
    assert row['indel_count'] == 1

    # Check last window (400k)
    row = result[result['window'] == 400_000].iloc[0]
    assert row['variant_count'] == 1


def test_compute_window_density_custom_window(sample_variant_df):
    result = compute_window_density(sample_variant_df, window_size=200_000)

    # Expected windows: 0, 200k, 400k
    assert len(result) == 3

    # Window 0 → positions 100k, 150k
    row0 = result[result['window'] == 0].iloc[0]
    assert row0['variant_count'] == 2

    # Window 200k → positions 250k, 300k
    row1 = result[result['window'] == 200_000].iloc[0]
    assert row1['variant_count'] == 2

    # Window 400k → position 450k
    row2 = result[result['window'] == 400_000].iloc[0]
    assert row2['variant_count'] == 1


# ─────────────────────────────────────────
# TEST: melt_population_af
# ─────────────────────────────────────────

def test_melt_population_af(sample_variant_df):
    result = melt_population_af(sample_variant_df)

    assert isinstance(result, pd.DataFrame)

    # One NaN already present → total = 25 - 1 = 24
    assert len(result) == 24

    assert list(result.columns) == ['pos', 'variant_type', 'population', 'allele_frequency']

    # Check population labels
    assert set(result['population'].unique()) == {
        'African', 'European', 'East Asian', 'American', 'South Asian'
    }

    # Ensure no NaNs
    assert result['allele_frequency'].isna().sum() == 0


def test_melt_population_af_drops_na(sample_variant_df):
    df_with_na = sample_variant_df.copy()
    df_with_na.loc[0, 'af_amr'] = None  # add one more NaN

    result = melt_population_af(df_with_na)

    # Originally 1 NaN → now 2 NaNs → 25 - 2 = 23
    assert len(result) == 23
