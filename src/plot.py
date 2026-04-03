# src/plot.py

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

def save_figure(fig, output_path, filename):
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path / filename, bbox_inches="tight")
    plt.close(fig)

def plot_variant_density(df, window_size=100_000, chrom=None, save=False, output_dir="outputs/figures"):
    data = df.copy()

    if chrom:
        data = data[data["chrom"] == chrom]

    data["window"] = (data["pos"] // window_size) * window_size
    density = data.groupby("window").size().reset_index(name="variant_count")

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.bar(density["window"] / 1e6, density["variant_count"], width=(window_size / 1e6) * 0.8)

    ax.set_xlabel("Genomic Window Start (Mb)")
    ax.set_ylabel("Variant Count")
    ax.set_title(f"Variant Density per {window_size // 1000}kb")

    if save:
        save_figure(fig, output_dir, "variant_density.png")
    else:
        plt.show()

def plot_af_violin(df, save=False, output_dir="outputs/figures"):
    pop_cols = ['af_afr', 'af_eur', 'af_eas', 'af_amr', 'af_sas']
    pop_labels = {
        'af_afr': 'African', 'af_eur': 'European', 
        'af_eas': 'East Asian', 'af_amr': 'American', 'af_sas': 'South Asian'
    }

    melted = df[['pos', 'variant_type'] + pop_cols].melt(
        id_vars=['pos', 'variant_type'],
        value_vars=pop_cols,
        var_name='population',
        value_name='AF'
    )
    melted['population'] = melted['population'].map(pop_labels)
    

    melted = melted.dropna(subset=['AF'])
    melted = melted[melted['AF'] > 0]

    fig, ax = plt.subplots(figsize=(10, 5))
    
    # inner="quartile" adds clean dashed lines inside the violin for median and interquartile ranges
    sns.violinplot(data=melted, x="population", y="AF", ax=ax, palette="muted", inner="quartile")

    ax.set_yscale('log')

    ax.set_xlabel("Population")
    ax.set_ylabel("Allele Frequency (AF) - Log Scale")
    ax.set_title("AF Distribution by Population (Log Scale)")

    if save:
        save_figure(fig, output_dir, "af_violin.png")
    else:
        plt.show()

def plot_variant_heatmap(df, region_size=1_000_000, save=False, output_dir="outputs/figures"):
    data = df.copy()
    # Create 1 Megabase (1Mb) regions to avoid plotting thousands of tiny windows
    data["region_Mb"] = (data["pos"] // region_size) * region_size / 1e6 
    
    heatmap_data = pd.pivot_table(
        data,
        index="region_Mb",
        columns="variant_type",
        aggfunc="size",
        fill_value=0
    )

    fig, ax = plt.subplots(figsize=(8, 10))
    sns.heatmap(heatmap_data, annot=True, fmt="d", ax=ax, cmap="YlGnBu")

    ax.set_xlabel("Variant Type")
    ax.set_ylabel("Genomic Region (Mb)")
    ax.set_title("Variant Type Distribution Across Regions")

    if save:
        save_figure(fig, output_dir, "variant_heatmap.png")
    else:
        plt.show()


def plot_af_scatter(df, pop1="EUR", pop2="AFR", save=False, output_dir="outputs/figures"):

    col1 = f"af_{pop1.lower()}" if not pop1.lower().startswith("af_") else pop1.lower()
    col2 = f"af_{pop2.lower()}" if not pop2.lower().startswith("af_") else pop2.lower()

    plot_df = df[[col1, col2]].dropna()

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(plot_df[col1], plot_df[col2], alpha=0.3, s=15, c="purple")

    ax.set_xlabel(f"Allele Frequency ({pop1})")
    ax.set_ylabel(f"Allele Frequency ({pop2})")
    ax.set_title(f"AF Comparison: {pop1} vs {pop2}")


    min_val = min(plot_df[col1].min(), plot_df[col2].min())
    max_val = max(plot_df[col1].max(), plot_df[col2].max())
    if pd.notna(min_val) and pd.notna(max_val):
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', alpha=0.5)

    if save:
        save_figure(fig, output_dir, f"af_scatter_{pop1}_{pop2}.png")
    else:
        plt.show()

import seaborn as sns
import matplotlib.pyplot as plt

def plot_af_histograms(df, save=False, output_dir="outputs/figures"):
    pop_cols = ['af_afr', 'af_eur', 'af_eas', 'af_amr', 'af_sas']
    pop_labels = {
        'af_afr': 'African', 'af_eur': 'European', 
        'af_eas': 'East Asian', 'af_amr': 'American', 'af_sas': 'South Asian'
    }

    melted = df[['pos', 'variant_type'] + pop_cols].melt(
        id_vars=['pos', 'variant_type'],
        value_vars=pop_cols,
        var_name='population',
        value_name='AF'
    )
    melted['population'] = melted['population'].map(pop_labels)
    
    melted = melted.dropna(subset=['AF'])
    melted = melted[melted['AF'] > 0]


    g = sns.displot(
        data=melted, 
        x="AF", 
        col="population",   
        col_wrap=3,         
        kind="hist",        
        log_scale=True,     
        bins=30,            
        height=4,           
        color="teal",
        facet_kws={'sharex': False, 'sharey': False} 
    )


    for ax in g.axes.flatten():
        if ax is not None:
            ax.tick_params(labelbottom=True)
            ax.set_xlabel("Allele Frequency (AF) - Log Scale")


    g.fig.suptitle("Allele Frequency Distribution by Population", y=1.02, fontsize=16)
    

    plt.tight_layout()

    if save:
        save_figure(g.fig, output_dir, "af_histograms_faceted.png")
    else:
        plt.show()
