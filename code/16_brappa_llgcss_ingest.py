#!/usr/bin/env python3
"""
16_brappa_llgcss_ingest.py
Ingest and harmonize the Brassica rapa (Lunar LEAF / LLGCSS) GCR radiation and floral scent datasets.
Maps B. rapa gene identifiers (BRA###### / Bra######) to Arabidopsis thaliana orthologs (TAIR10 AGI),
annotates curated floral volatile pathways, and formats query signatures for cross-species decoding.

Author: Antigravity Assistant & Dr. Richard Barker
Date: November 2026
"""

import os
import sys
from pathlib import Path
import pandas as pd
import numpy as np

# Base paths
WORKSPACE_DIR = Path(__file__).resolve().parents[1]
LLGCSS_DIR = Path("/Users/drb_laptop/Documents/B_rappa_LLGCSS")
DATA_OUT = WORKSPACE_DIR / "data"
DATA_OUT.mkdir(exist_ok=True)

print("=" * 70)
print("16_brappa_llgcss_ingest: Starting B. rapa LLGCSS data ingestion")
print("=" * 70)

# 1. Load Orthology Map (18,770 RBH pairs)
ortho_file = LLGCSS_DIR / "annotation" / "brapa_to_arabidopsis_orthologs.tsv"
print(f"Loading orthology map from: {ortho_file}")
df_ortho = pd.read_csv(ortho_file, sep="\t")
# Columns: brapa_gene, arabidopsis_agi, bitscore, identity, evalue, qcov, scov
df_ortho["BRA_id"] = df_ortho["brapa_gene"].str.strip().str.upper()
df_ortho["bra_id"] = df_ortho["brapa_gene"].str.strip().str.lower().str.capitalize()
df_ortho["arabidopsis_agi"] = df_ortho["arabidopsis_agi"].str.strip().str.upper()
print(f"Loaded {len(df_ortho):,} 1:1 ortholog pairs between B. rapa and A. thaliana.")

# 2. Load Curated Scent Gene Set (363 genes)
scent_geneset_file = LLGCSS_DIR / "annotation" / "scent_geneset.tsv"
print(f"Loading curated scent gene set from: {scent_geneset_file}")
df_scent_set = pd.read_csv(scent_geneset_file, sep="\t")
# Columns: bra_id, BRA_id, in_scent_data, in_radiation_data, tier, routes, evidence
df_scent_set["BRA_id"] = df_scent_set["BRA_id"].str.strip().str.upper()
df_scent_set = df_scent_set[["BRA_id", "tier", "routes", "evidence"]].drop_duplicates("BRA_id")
print(f"Loaded {len(df_scent_set)} curated scent genes.")

# 3. Load Radiation Model Results (DESeq2 contrasts)
rad_results_file = LLGCSS_DIR / "Radiation model" / "Results_LFC_Pval_DESeq2.csv"
print(f"Loading radiation contrasts from: {rad_results_file}")
df_rad = pd.read_csv(rad_results_file)
# Normalize columns and gene id
df_rad["BRA_id"] = df_rad["ensembl_ID"].astype(str).str.strip().str.upper()

# 4. Load Interaction Terms Model
rad_inter_file = LLGCSS_DIR / "Radiation model" / "Radi_model_with_interactions_deg_values_DESeq2.csv"
print(f"Loading interaction terms from: {rad_inter_file}")
df_inter = pd.read_csv(rad_inter_file)
df_inter["BRA_id"] = df_inter["ensembl_ID"].astype(str).str.strip().str.upper()

# Select essential columns from radiation model
rad_cols = {
    "BRA_id": "BRA_id",
    "baseMean": "baseMean",
    "WT_&_40_CGY-WT_&_0_CGY_log2FC": "WT_40cGy_vs_0cGy_log2FC",
    "WT_&_40_CGY-WT_&_0_CGY_adjPval": "WT_40cGy_vs_0cGy_padj",
    "ANTHO_LESS_&_40_CGY-ANTHO_LESS_&_0_CGY_log2FC": "antho_less_40cGy_vs_0cGy_log2FC",
    "ANTHO_LESS_&_40_CGY-ANTHO_LESS_&_0_CGY_adjPval": "antho_less_40cGy_vs_0cGy_padj",
    "ANTHO_LESS_&_0_CGY-WT_&_0_CGY_log2FC": "antho_less_vs_WT_0cGy_log2FC",
    "ANTHO_LESS_&_0_CGY-WT_&_0_CGY_adjPval": "antho_less_vs_WT_0cGy_padj",
    "ANTHO_LESS_&_40_CGY-WT_&_40_CGY_log2FC": "antho_less_vs_WT_40cGy_log2FC",
    "ANTHO_LESS_&_40_CGY-WT_&_40_CGY_adjPval": "antho_less_vs_WT_40cGy_padj"
}
rad_extracted = df_rad[[c for c in rad_cols.keys() if c in df_rad.columns]].rename(columns=rad_cols)

# Select interaction columns
inter_cols = {
    "BRA_id": "BRA_id",
    "I:condition_40_CGY.genotype_ANTHO_LESS_log2FC": "inter_dose_x_genotype_log2FC",
    "I:condition_40_CGY.genotype_ANTHO_LESS_adjPval": "inter_dose_x_genotype_padj",
    "I:condition_40_CGY.preservative_RL_log2FC": "inter_dose_x_preservative_log2FC",
    "I:condition_40_CGY.preservative_RL_adjPval": "inter_dose_x_preservative_padj"
}
inter_extracted = df_inter[[c for c in inter_cols.keys() if c in df_inter.columns]].rename(columns=inter_cols)

# Merge radiation contrasts
rad_merged = pd.merge(rad_extracted, inter_extracted, on="BRA_id", how="outer")
print(f"Radiation dataset contains {len(rad_merged):,} genes.")

# 5. Load Floral Scent High vs Low Results
scent_de_file = LLGCSS_DIR / "NewTest" / "differential_expression_GLbulkRNAseq (1).csv"
print(f"Loading scent DE results from: {scent_de_file}")
df_scent_de = pd.read_csv(scent_de_file)
df_scent_de["BRA_id"] = df_scent_de["gene_id"].astype(str).str.strip().str.upper()

scent_extracted = df_scent_de[[
    "BRA_id",
    "Log2fc_(High)v(Low)",
    "Adj.p.value_(High)v(Low)",
    "All.mean"
]].rename(columns={
    "Log2fc_(High)v(Low)": "scent_High_vs_Low_log2FC",
    "Adj.p.value_(High)v(Low)": "scent_High_vs_Low_padj",
    "All.mean": "scent_baseMean"
})
print(f"Scent dataset contains {len(scent_extracted):,} genes.")

# 6. Master Integration: Orthology + Scent Annotation + GCR Radiation + Scent DE
master = pd.merge(rad_merged, df_ortho, on="BRA_id", how="left")
master = pd.merge(master, df_scent_set, on="BRA_id", how="left")
master = pd.merge(master, scent_extracted, on="BRA_id", how="left")

# Clean metadata tags
master["has_arabidopsis_ortholog"] = ~master["arabidopsis_agi"].isna()
master["is_scent_curated"] = ~master["tier"].isna()

print(f"\nMaster table summary:")
print(f"  Total genes: {len(master):,}")
print(f"  Genes with Arabidopsis TAIR10 ortholog: {master['has_arabidopsis_ortholog'].sum():,}")
print(f"  Curated floral scent genes present: {master['is_scent_curated'].sum():,}")

# Export master dataset
out_master_path = DATA_OUT / "brappa_master_orthology_contrasts.csv"
master.to_csv(out_master_path, index=False)
print(f"\n[Saved] Master table -> {out_master_path}")

# Export specific query signatures for DSRS and Redox CVAE
# 1. WT 40 cGy vs 0 cGy (Arabidopsis AGI keyed)
sig_wt = master.dropna(subset=["arabidopsis_agi", "WT_40cGy_vs_0cGy_log2FC"])[
    ["arabidopsis_agi", "WT_40cGy_vs_0cGy_log2FC", "WT_40cGy_vs_0cGy_padj"]
].drop_duplicates("arabidopsis_agi").rename(columns={
    "arabidopsis_agi": "gene_id",
    "WT_40cGy_vs_0cGy_log2FC": "log2FC",
    "WT_40cGy_vs_0cGy_padj": "padj"
})
sig_wt_path = DATA_OUT / "brappa_sig_wt_gcr_40cGy.csv"
sig_wt.to_csv(sig_wt_path, index=False)
print(f"[Saved] WT 40 cGy GCR signature ({len(sig_wt):,} genes) -> {sig_wt_path}")

# 2. Anthocyaninless 40 cGy vs 0 cGy
sig_anl = master.dropna(subset=["arabidopsis_agi", "antho_less_40cGy_vs_0cGy_log2FC"])[
    ["arabidopsis_agi", "antho_less_40cGy_vs_0cGy_log2FC", "antho_less_40cGy_vs_0cGy_padj"]
].drop_duplicates("arabidopsis_agi").rename(columns={
    "arabidopsis_agi": "gene_id",
    "antho_less_40cGy_vs_0cGy_log2FC": "log2FC",
    "antho_less_40cGy_vs_0cGy_padj": "padj"
})
sig_anl_path = DATA_OUT / "brappa_sig_antho_less_gcr_40cGy.csv"
sig_anl.to_csv(sig_anl_path, index=False)
print(f"[Saved] Antho_less 40 cGy GCR signature ({len(sig_anl):,} genes) -> {sig_anl_path}")

# 3. Interaction Condition 40 cGy x Genotype Antho_less
sig_inter = master.dropna(subset=["arabidopsis_agi", "inter_dose_x_genotype_log2FC"])[
    ["arabidopsis_agi", "inter_dose_x_genotype_log2FC", "inter_dose_x_genotype_padj"]
].drop_duplicates("arabidopsis_agi").rename(columns={
    "arabidopsis_agi": "gene_id",
    "inter_dose_x_genotype_log2FC": "log2FC",
    "inter_dose_x_genotype_padj": "padj"
})
sig_inter_path = DATA_OUT / "brappa_sig_inter_dose_genotype.csv"
sig_inter.to_csv(sig_inter_path, index=False)
print(f"[Saved] Interaction dose x genotype signature ({len(sig_inter):,} genes) -> {sig_inter_path}")

# 4. Scent High vs Low
sig_scent = master.dropna(subset=["arabidopsis_agi", "scent_High_vs_Low_log2FC"])[
    ["arabidopsis_agi", "scent_High_vs_Low_log2FC", "scent_High_vs_Low_padj"]
].drop_duplicates("arabidopsis_agi").rename(columns={
    "arabidopsis_agi": "gene_id",
    "scent_High_vs_Low_log2FC": "log2FC",
    "scent_High_vs_Low_padj": "padj"
})
sig_scent_path = DATA_OUT / "brappa_sig_scent_high_vs_low.csv"
sig_scent.to_csv(sig_scent_path, index=False)
print(f"[Saved] Scent High vs Low signature ({len(sig_scent):,} genes) -> {sig_scent_path}")

print("=" * 70)
print("16_brappa_llgcss_ingest: Completed successfully!")
print("=" * 70)
