#!/usr/bin/env python3
"""
Visualization functions for PRISTINE analysis results.

This module contains plotting functions for generating heatmaps and
distribution plots of SNP informativeness and density.
"""

from __future__ import annotations
import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


def plot_informativeness_heatmap(csv_file, output_file):
    """
    Generate a heatmap showing informative SNP proportions per locus and species.

    :param csv_file: Path to the SNP summary CSV file
    :param output_file: Path where the heatmap PNG will be saved
    :return: Path to the saved heatmap file
    """
    df = pd.read_csv(csv_file)

    prop_cols = [col for col in df.columns if col.startswith("Prop_") and col != "Avg_Prop_Informative_SNPs"]
    if not prop_cols:
        print("No proportion columns found in the dataset.")
        return None

    heatmap_data = df.set_index("Locus")[prop_cols]

    if "Avg_Prop_Informative_SNPs" in df.columns:
        sorted_loci = df.sort_values("Avg_Prop_Informative_SNPs", ascending=False)["Locus"]
        heatmap_data = heatmap_data.loc[sorted_loci]

    plt.figure(figsize=(12, max(6, 0.3 * len(heatmap_data))))
    sns.heatmap(heatmap_data, annot=False, cmap="YlOrRd", cbar_kws={'label': 'Proportion of Informative SNPs'})
    plt.title("Informative SNP Proportions per Locus and Non-Target Species")
    plt.xlabel("Non-Target Species")
    plt.ylabel("Locus")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.close()

    return output_file


def plot_snp_density_distribution(csv_file, output_dir, top_n=5):
    """
    Generate SNP position distribution plots for top-ranked loci.

    :param csv_file: Path to the SNP summary CSV file
    :param output_dir: Directory where SNP density plots will be saved
    :param top_n: Number of top-ranked loci to plot (default: 5)
    :return: Path to the directory containing the plots
    """
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(csv_file)

    if "Avg_Prop_Informative_SNPs" not in df.columns:
        print("Missing ranking column 'Avg_Prop_Informative_SNPs'.")
        return None

    top_loci = df.sort_values("Avg_Prop_Informative_SNPs", ascending=False).head(top_n)

    prop_cols = [col for col in df.columns if col.startswith("SNP_Pos_")]

    for _, row in top_loci.iterrows():
        locus_name = row["Locus"]
        plt.figure(figsize=(12, 3))

        for col in prop_cols:
            species = col.replace("SNP_Pos_", "")
            if pd.isna(row[col]) or not row[col].strip():
                continue
            try:
                positions = list(map(int, row[col].split(',')))
            except ValueError:
                continue
            sns.histplot(positions, bins=50, kde=False, label=species, element="step", fill=False)

        plt.title(f"SNP Position Distribution - {locus_name}")
        plt.xlabel("Alignment Position")
        plt.ylabel("SNP Count")
        plt.legend(title="Non-Target Species", loc="upper right", fontsize="small")
        plt.tight_layout()
        plot_path = os.path.join(output_dir, f"snp_density_{locus_name.replace('.', '_')}.png")
        plt.savefig(plot_path)
        plt.close()

    print(f"SNP density plots saved to: {output_dir}")
    return output_dir
