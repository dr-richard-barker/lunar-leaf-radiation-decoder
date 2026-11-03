#!/usr/bin/env python3
"""
19_brappa_rri_and_scent_tradeoff.py
Phase 3 — Compute the Radiation Resilience Index (RRI) across the 31 B. rapa
libraries and quantify the metabolic trade-off between DNA repair/antioxidant defense
and floral volatile emission.

Evaluates:
1. Sample-level RRI (latent stability, pathway balance, and housekeeping module preservation)
2. Genotype comparison: WT vs anthocyaninless mutant resilience under 40 cGy GCR
3. Preservative effect: DRS vs RNAlater stability
4. Scent pathway routing: Benzenoid/phenylpropanoid, terpenoid, and carotenoid trade-offs
5. DSB repair core activation (BRCA1, PARP2, RAD51, SMR7) vs volatile synthases

Author: Antigravity Assistant & Dr. Richard Barker
Date: November 2026
"""

import os
import sys
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.stats import entropy, spearmanr, pearsonr

WORKSPACE_DIR = Path(__file__).resolve().parents[1]
DATA_DIR = WORKSPACE_DIR / "data"
RESULTS_DIR = WORKSPACE_DIR / "results"
RESULTS_DIR.mkdir(exist_ok=True)

print("=" * 70)
print("19_brappa_rri_and_scent_tradeoff: Computing RRI & Metabolic Trade-offs")
print("=" * 70)

# 1. Load Expression Data & Metadata from B. rapa
deseq_file = Path("/Users/drb_laptop/Documents/B_rappa_LLGCSS/Radiation model/Results_LFC_Pval_DESeq2.csv")
meta_file = Path("/Users/drb_laptop/Documents/B_rappa_LLGCSS/counts and factors/LEAF_Metadata.csv")
master_file = DATA_DIR / "brappa_master_orthology_contrasts.csv"

print(f"Loading DESeq2 expression matrix from: {deseq_file}")
df_deseq = pd.read_csv(deseq_file)
df_meta = pd.read_csv(meta_file)
df_master = pd.read_csv(master_file)

# Identify sample columns
sample_cols = [c for c in df_deseq.columns if "GLbulkRNAseq" in c]
print(f"Identified {len(sample_cols)} RNA-seq sample libraries in DESeq2 matrix.")

# Convert sample columns to float and fill NaN with 0.0
for c in sample_cols:
    df_deseq[c] = pd.to_numeric(df_deseq[c], errors="coerce").fillna(0.0)

# Clean gene IDs
df_deseq["BRA_id"] = df_deseq["ensembl_ID"].astype(str).str.strip().str.upper()

# Merge with ortholog and scent annotations
expr_annot = pd.merge(
    df_deseq[["BRA_id"] + sample_cols],
    df_master[["BRA_id", "arabidopsis_agi", "tier", "routes", "is_scent_curated", "WT_40cGy_vs_0cGy_log2FC", "antho_less_40cGy_vs_0cGy_log2FC"]],
    on="BRA_id",
    how="left"
)

# 2. Extract Subsets for RRI computation
# Filter for genes with detectable expression across samples
row_sums = expr_annot[sample_cols].sum(axis=1)
expr_expressed = expr_annot[row_sums > 0].copy()

# A. Ortholog space for latent distance
ortho_expr = expr_expressed.dropna(subset=["arabidopsis_agi"]).drop_duplicates("arabidopsis_agi")
X = ortho_expr[sample_cols].values.T # (31 samples x N genes)

# Parse metadata per sample
meta_records = []
for s in sample_cols:
    row = df_meta[df_meta["sample_id"] == s]
    if not row.empty:
        cond = row["condition"].values[0]
        pres = row["preservative"].values[0]
        geno = row["genotype"].values[0]
    else:
        # Infer from name
        cond = "40_cGy" if "40" in s else "0_cGy"
        pres = "RL" if "R" in s else "DRS"
        geno = "antho_less" if s.startswith("AL") else "WT"
    is_ctrl = (cond == "0_cGy") and (geno == "WT")
    meta_records.append({
        "sample_id": s,
        "condition": cond,
        "preservative": pres,
        "genotype": geno,
        "is_control": is_ctrl
    })
meta_df = pd.DataFrame(meta_records)

# 3. Compute RRI Components
# Component 1: Latent Distance RRI
ctrl_mask = meta_df["is_control"].values
if np.sum(ctrl_mask) == 0:
    ctrl_mask = (meta_df["condition"] == "0_cGy").values

ctrl_X = X[ctrl_mask]
ctrl_centroid = np.mean(ctrl_X, axis=0, keepdims=True)

# Pairwise distances among controls for sigma
ctrl_dists = cdist(ctrl_X, ctrl_X, metric="euclidean")
sigma = np.median(ctrl_dists[ctrl_dists > 0]) or 1.0

# Distance of each sample to control centroid
sample_dists = cdist(X, ctrl_centroid, metric="euclidean").flatten()
rri_latent = np.exp(- (sample_dists ** 2) / (2 * (sigma ** 2)))
rri_latent = np.clip(rri_latent, 0.0, 1.0)

# Component 2: Pathway Balance RRI
# Pathway sets: DSB Repair, ROS Scavenging, Anthocyanin/Phenylpropanoid, Scent/Terpenoid
pathway_genes = {
    "dsb_repair": ["AT4G21070", "AT4G02390", "AT5G20850", "AT3G27630", "AT3G02820", "AT5G41650"],
    "ros_defense": ["AT1G07890", "AT4G35090", "AT1G20630", "AT1G08830", "AT4G25100", "AT3G22370", "AT3G24170"],
    "phenylpropanoid": ["AT5G13930", "AT3G53260", "AT5G42800", "AT3G51240"]
}

pathway_scores_per_sample = []
for s in sample_cols:
    scores = []
    # DSB
    sub_dsb = expr_expressed[expr_expressed["arabidopsis_agi"].isin(pathway_genes["dsb_repair"])]
    scores.append(sub_dsb[s].mean() if not sub_dsb.empty else 1.0)
    # ROS
    sub_ros = expr_expressed[expr_expressed["arabidopsis_agi"].isin(pathway_genes["ros_defense"])]
    scores.append(sub_ros[s].mean() if not sub_ros.empty else 1.0)
    # Scent Tier 1
    sub_scent = expr_expressed[expr_expressed["tier"] == 1]
    scores.append(sub_scent[s].mean() if not sub_scent.empty else 1.0)
    
    # Normalize to probabilities for Shannon evenness
    p_vec = np.array(scores)
    p_norm = p_vec / np.sum(p_vec) if np.sum(p_vec) > 0 else np.ones(len(scores)) / len(scores)
    h = entropy(p_norm)
    h_max = np.log(len(scores))
    evenness = h / h_max if h_max > 0 else 1.0
    pathway_scores_per_sample.append(evenness)

rri_pathway = np.clip(np.array(pathway_scores_per_sample), 0.0, 1.0)

# Component 3: Housekeeping Module Preservation RRI
# Use top 20% lowest-variance genes across controls as proxy for sustained housekeeping module
ctrl_var = np.var(ctrl_X, axis=0)
hk_idx = np.argsort(ctrl_var)[:int(0.2 * len(ctrl_var))]
ctrl_hk_mean = np.mean(ctrl_X[:, hk_idx], axis=0)

rri_module = []
for i in range(len(sample_cols)):
    r, _ = pearsonr(X[i, hk_idx], ctrl_hk_mean)
    rri_module.append(max(0.0, r))
rri_module = np.array(rri_module)

# Composite RRI: 0.50 latent + 0.25 pathway + 0.25 module
rri_total = 0.50 * rri_latent + 0.25 * rri_pathway + 0.25 * rri_module

meta_df["RRI"] = rri_total
meta_df["RRI_latent"] = rri_latent
meta_df["RRI_pathway"] = rri_pathway
meta_df["RRI_module"] = rri_module

print(f"\nRRI Summary Statistics:")
rri_summary = meta_df.groupby(["genotype", "condition"])[["RRI", "RRI_latent", "RRI_pathway", "RRI_module"]].mean()
print(rri_summary)

# 4. Scent vs Radiation Trade-off Analysis
# Group scent genes by route and test expression difference across WT vs AnthoLess under 40 cGy
scent_df = expr_annot[expr_annot["is_scent_curated"] == True].copy()
tradeoff_results = []

wt_0_cols = [s for s in sample_cols if "W0" in s]
wt_40_cols = [s for s in sample_cols if "W40" in s]
al_0_cols = [s for s in sample_cols if "AL0" in s]
al_40_cols = [s for s in sample_cols if "AL40" in s]

for route, grp in scent_df.groupby("routes"):
    wt_0_vals = grp[wt_0_cols].values.flatten()
    wt_40_vals = grp[wt_40_cols].values.flatten()
    al_0_vals = grp[al_0_cols].values.flatten()
    al_40_vals = grp[al_40_cols].values.flatten()
    
    wt_fold = np.mean(wt_40_vals) - np.mean(wt_0_vals)
    al_fold = np.mean(al_40_vals) - np.mean(al_0_vals)
    
    tradeoff_results.append({
        "route": route,
        "n_genes": len(grp),
        "wt_mean_delta": wt_fold,
        "al_mean_delta": al_fold,
        "delta_gap": al_fold - wt_fold
    })

df_tradeoff = pd.DataFrame(tradeoff_results).sort_values("n_genes", ascending=False)
print(f"\nFloral Scent Route Metabolic Shifts (40 cGy vs 0 cGy):")
print(df_tradeoff)

# Save RRI and Trade-off Tables
out_rri_path = RESULTS_DIR / "brappa_rri_sample_scores.csv"
out_tradeoff_path = RESULTS_DIR / "brappa_scent_tradeoff_summary.csv"

meta_df.to_csv(out_rri_path, index=False)
df_tradeoff.to_csv(out_tradeoff_path, index=False)

print(f"\n[Saved] Sample RRI Scores -> {out_rri_path}")
print(f"[Saved] Scent Tradeoff Summary -> {out_tradeoff_path}")
print("=" * 70)
print("19_brappa_rri_and_scent_tradeoff: Completed successfully!")
print("=" * 70)
