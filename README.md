# Genomic Variant Landscape & Population Comparison Dashboard

## Overview
This project analyzes genomic variant patterns using public human genetic variation data. The goal is to explore how variant density and allele frequencies vary across genomic regions and populations. The analysis focuses on handling large genomic datasets and producing interpretable visual summaries.

## Data Source
The dataset is obtained from the Genome Aggregation Database (gnomAD).

Dataset used:
- `gnomad.genomes.v4.1.sites.chrY.vcf.bgz`
- Contains variant sites on the human Y chromosome.

The file includes:
- genomic coordinates (chromosome and position)
- reference and alternate alleles
- variant filtering status
- allele counts and allele frequencies
- population-specific frequency summaries

## Data Processing

### 1. Data Ingestion
The compressed VCF file (`.vcf.bgz`) is read directly in Python and parsed into a dataframe.

### 2. Parsing Variant Information
The `INFO` column contains multiple annotations (such as allele counts and frequencies). These are extracted and converted into separate columns for easier analysis.

### 3. Filtering
Variants can be filtered based on:
- filter status (`PASS`)
- selected genomic regions
- allele frequency thresholds

## Analysis
The project focuses on summarizing genomic variation across regions and populations.

Examples include:
- Variant density across genomic windows along the Y chromosome
- Allele frequency distributions
- Population-specific allele frequency comparisons

These analyses help identify genomic regions with unusually high or low variant density.

## Data Representation Comparison
Two representations of the data are compared:

1. **Raw compressed VCF**
   - Genomic format optimized for storage and bioinformatics tools.

2. **Parsed tabular dataset**
   - INFO fields converted into structured columns.
   - Easier to analyze and visualize in Python.

This comparison highlights trade-offs in usability and processing efficiency.
###
