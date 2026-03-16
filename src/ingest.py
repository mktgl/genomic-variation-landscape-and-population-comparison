# src/ingest.py
# ─────────────────────────────────────────────────────────────
# PURPOSE: Functions to read and parse VCF files into DataFrames
# Used by: notebooks/01_data_ingestion.ipynb
# ─────────────────────────────────────────────────────────────

import cyvcf2
import pandas as pd


def parse_vcf_to_dataframe(vcf_path: str, max_rows: int = None) -> pd.DataFrame:
    """
    Reads a .vcf or .vcf.bgz file and returns a pandas DataFrame.

    Parameters:
        vcf_path : path to the VCF file
        max_rows : if set, stops after this many variants (useful for testing)

    Returns:
        DataFrame with columns: chrom, pos, ref, alt, qual, filter_status,
                                 af_global, af_afr, af_eur, af_eas, af_amr, af_sas,
                                 ac, an, variant_type
    """
    vcf = cyvcf2.VCF(vcf_path)
    rows = []

    for i, variant in enumerate(vcf):
        if max_rows and i >= max_rows:
            break

        rows.append({
            'chrom'        : variant.CHROM,
            'pos'          : variant.POS,
            'ref'          : variant.REF,
            'alt'          : str(variant.ALT[0]) if variant.ALT else None,
            'qual'         : variant.QUAL,
            'filter_status': variant.FILTER if variant.FILTER else 'PASS',
            'af_global'    : variant.INFO.get('AF'),
            'af_afr'       : variant.INFO.get('AF_afr'),
            'af_eur'       : variant.INFO.get('AF_nfe'),
            'af_eas'       : variant.INFO.get('AF_eas'),
            'af_amr'       : variant.INFO.get('AF_amr'),
            'af_sas'       : variant.INFO.get('AF_sas'),
            'ac'           : variant.INFO.get('AC'),
            'an'           : variant.INFO.get('AN'),
            'variant_type' : 'SNP' if (len(variant.REF) == 1 and
                                        variant.ALT and
                                        len(variant.ALT[0]) == 1) else 'INDEL'
        })

    vcf.close()
    return pd.DataFrame(rows)
