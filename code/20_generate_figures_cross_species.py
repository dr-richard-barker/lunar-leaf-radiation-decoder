#!/usr/bin/env python3
"""
20_generate_figures_cross_species.py
Phase 4 — Generate publication-quality figures (Figures 1–6) for the
Cross-Species Radiation & Redox Decoding manuscript.

Figures produced:
- Fig 1: Conceptual architecture & the classification vs mechanism continuum
- Fig 2: DSRS Stress Recognition & Cross-Species Spaceflight Concordance
- Fig 3: Redox CVAE Decoding & The "Antioxidant Shield Deficit" (WT vs antho_less)
- Fig 4: Radiation Resilience Index (RRI) Degradation & Housekeeping Collapse
- Fig 5: Metabolic Trade-off: Floral Scent Volatiles vs DNA Double-Strand Break Repair
- Fig 6: Spatial Anatomical Vulnerability Across Leaf, Root, and Floral Tissues

Author: Antigravity Assistant & Dr. Richard Barker
Date: November 2026
"""

import os
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns

WORKSPACE_DIR = Path(__file__).resolve().parents[1]
RESULTS_DIR = WORKSPACE_DIR / "results"
FIG_DIR = RESULTS_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True, parents=True)

# Styling configuration
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
    "font.size": 10,
    "axes.labelsize": 11,
    "axes.titlesize": 12,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.titlesize": 13,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight"
})

print("=" * 70)
print("20_generate_figures_cross_species: Generating Publication Figures")
print("=" * 70)

# Load results data
dsrs_ranks = pd.read_csv(RESULTS_DIR / "dsrs_brappa_stressor_rankings.csv")
dsrs_fam = pd.read_csv(RESULTS_DIR / "dsrs_brappa_family_rankings.csv")
redox_preds = pd.read_csv(RESULTS_DIR / "redox_decoder_brappa_predictions.csv")
redox_spat = pd.read_csv(RESULTS_DIR / "redox_decoder_organ_cell_scores.csv")
rri_scores = pd.read_csv(RESULTS_DIR / "brappa_rri_sample_scores.csv")
scent_tradeoff = pd.read_csv(RESULTS_DIR / "brappa_scent_tradeoff_summary.csv")

# -------------------------------------------------------------------------
# Figure 1: Conceptual Architecture & The Classification Continuum
# -------------------------------------------------------------------------
print("Generating Figure 1: Conceptual Architecture...")
fig, ax = plt.subplots(figsize=(10, 6))
ax.axis("off")

# Draw hierarchical boxes
boxes = [
    {"title": "1. Input Cohort (Lunar LEAF)", "text": "Brassica rapa (Wisconsin Fast Plant)\n• 40 cGy GCR radiation (WT vs antho_less)\n• Floral scent lines (High vs Low)", "xy": (0.05, 0.65), "w": 0.26, "h": 0.28, "color": "#E8F0FE", "edge": "#1A73E8"},
    {"title": "2. Orthology Bridge", "text": "18,770 1:1 Reciprocal Best Hits\n(DIAMOND proteome alignment)\nBra###### / BRA###### → TAIR10 AGI", "xy": (0.37, 0.65), "w": 0.26, "h": 0.28, "color": "#E6F4EA", "edge": "#1E8E3E"},
    {"title": "3. The Diagnostic Decoder Stack", "text": "Tier 1: DSRS Stress Recognition (27 stressors)\nTier 2: Redox CVAE Disentanglement (15 ROS classes)\nTier 3: GP-AE Latent & Resilience Index (RRI)\nTier 4: Floral Scent vs Repair Trade-off", "xy": (0.69, 0.65), "w": 0.28, "h": 0.28, "color": "#FEF7E0", "edge": "#F9AB00"},
    {"title": "Biological Hypothesis & Mechanistic Synthesis", "text": "Is it 'just a classification exercise'?\n• Operational Diagnostic: Fast, automated stress-fingerprint recognition & duration estimation.\n• Mechanistic Translation: Validating that Arabidopsis latent manifolds generalize to crops;\n  proving that anthocyanin antioxidant loss triggers superoxide collapse and impairs floral scent.", "xy": (0.15, 0.12), "w": 0.70, "h": 0.38, "color": "#FCE8E6", "edge": "#D93025"}
]

for b in boxes:
    rect = patches.FancyBboxPatch(b["xy"], b["w"], b["h"], boxstyle="round,pad=0.02",
                                  facecolor=b["color"], edgecolor=b["edge"], linewidth=2)
    ax.add_patch(rect)
    ax.text(b["xy"][0] + b["w"]/2, b["xy"][1] + b["h"] - 0.05, b["title"],
            ha="center", va="center", weight="bold", fontsize=11, color="#202124")
    ax.text(b["xy"][0] + b["w"]/2, b["xy"][1] + b["h"]/2 - 0.03, b["text"],
            ha="center", va="center", fontsize=9, color="#3C4043")

# Draw connecting arrows
arrow_style = dict(facecolor="#5F6368", edgecolor="#5F6368", width=1.5, headwidth=8)
ax.annotate("", xy=(0.37, 0.79), xytext=(0.31, 0.79), arrowprops=arrow_style)
ax.annotate("", xy=(0.69, 0.79), xytext=(0.63, 0.79), arrowprops=arrow_style)
ax.annotate("", xy=(0.50, 0.50), xytext=(0.50, 0.65), arrowprops=arrow_style)

fig.suptitle("Figure 1: Cross-Species Framework for Spaceflight Stress & Redox Decoding", fontsize=13, weight="bold", y=0.98)
plt.tight_layout()
fig.savefig(FIG_DIR / "fig1_conceptual_architecture.png")
fig.savefig(FIG_DIR / "fig1_conceptual_architecture.pdf")
plt.close(fig)

# -------------------------------------------------------------------------
# Figure 2: DSRS Stress Recognition Across Spaceflight Contrasts
# -------------------------------------------------------------------------
print("Generating Figure 2: DSRS Stress Recognition...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.2, 1]})

# Panel A: Top Stressors for WT vs AnthoLess Interaction
inter_ranks = dsrs_ranks[dsrs_ranks["query"] == "Brapa_Interaction_Dose_x_Genotype"].head(10).copy()
colors = ["#D93025" if "rad" in s or "gcr" in s else "#1A73E8" for s in inter_ranks["stressor_class"]]
ax1.barh(inter_ranks["stressor"], inter_ranks["spearman_rho"], color=colors, edgecolor="black", alpha=0.85)
ax1.set_xlabel("Spearman Rank Correlation (ρ with Reference Stressor)")
ax1.set_title("A. DSRS Top Matches for Radiation × Genotype Interaction", weight="bold")
ax1.invert_yaxis()
ax1.grid(axis="x", linestyle="--", alpha=0.5)

# Panel B: Stressor Family Comparison across Queries
fam_pivot = dsrs_fam.pivot(index="stressor_class", columns="query", values="spearman_rho")
col_map = {
    "Brapa_WT_GCR_40cGy": "WT GCR 40cGy",
    "Brapa_AnthoLess_GCR_40cGy": "AnthoLess GCR 40cGy",
    "Brapa_Interaction_Dose_x_Genotype": "Interaction (Dose × Genotype)",
    "Brapa_Scent_High_vs_Low": "Floral Scent (High vs Low)"
}
fam_pivot = fam_pivot.rename(columns=col_map)
fam_pivot = fam_pivot.sort_values(by="Interaction (Dose × Genotype)", ascending=False)

sns.heatmap(fam_pivot, annot=True, fmt=".2f", cmap="vlag", center=0, ax=ax2, cbar_kws={"label": "Mean Spearman ρ"})
ax2.set_title("B. Stressor Family Fingerprint Concordance", weight="bold")
ax2.set_ylabel("Stressor Class / Family")

fig.suptitle("Figure 2: DSRS Deep-Space Stress Recognition Decodes GCR Radiation in Brassica rapa", fontsize=13, weight="bold", y=0.98)
plt.tight_layout()
fig.savefig(FIG_DIR / "fig2_dsrs_stress_recognition.png")
fig.savefig(FIG_DIR / "fig2_dsrs_stress_recognition.pdf")
plt.close(fig)

# -------------------------------------------------------------------------
# Figure 3: Redox CVAE Decoding & The Antioxidant Shield Deficit
# -------------------------------------------------------------------------
print("Generating Figure 3: Redox CVAE Decoding...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

# Panel A: ROS Stimulus Class Probabilities
stim_cols = [c for c in redox_preds.columns if c.startswith("prob_")]
stim_labels = [c.replace("prob_", "").replace("_", " ").title() for c in stim_cols]
pred_plot = redox_preds.set_index("contrast")[stim_cols].T
pred_plot.index = stim_labels

x = np.arange(len(stim_labels))
width = 0.2
ax1.bar(x - 1.5*width, pred_plot["WT_GCR_40cGy"] * 100, width, label="WT GCR 40cGy", color="#1E8E3E", edgecolor="black")
ax1.bar(x - 0.5*width, pred_plot["AnthoLess_GCR_40cGy"] * 100, width, label="AnthoLess GCR 40cGy", color="#D93025", edgecolor="black")
ax1.bar(x + 0.5*width, pred_plot["Interaction_Dose_x_Genotype"] * 100, width, label="Interaction Term", color="#F9AB00", edgecolor="black")
ax1.bar(x + 1.5*width, pred_plot["Scent_High_vs_Low"] * 100, width, label="Scent High vs Low", color="#9334E6", edgecolor="black")

ax1.set_xticks(x)
ax1.set_xticklabels(stim_labels, rotation=25, ha="right")
ax1.set_ylabel("Predicted Probability (%)")
ax1.set_title("A. Decoded Reactive Oxygen Species (ROS) Profiles", weight="bold")
ax1.legend(loc="upper right")
ax1.grid(axis="y", linestyle="--", alpha=0.5)

# Panel B: Inferred Kinetics Timing
contrasts_disp = ["WT GCR 40cGy", "AnthoLess GCR 40cGy", "Interaction Term", "Scent High vs Low"]
hours = redox_preds["inferred_hours"].values
colors_time = ["#1E8E3E", "#D93025", "#F9AB00", "#9334E6"]

bars = ax2.bar(contrasts_disp, hours, color=colors_time, edgecolor="black", width=0.55)
ax2.axhspan(0, 1.0, color="#E8F0FE", alpha=0.4, label="Immediate (<1h)")
ax2.axhspan(1.0, 4.0, color="#FEF7E0", alpha=0.4, label="Early Acute (1-4h)")
ax2.axhspan(4.0, 12.0, color="#FCE8E6", alpha=0.4, label="Mid Response (4-12h)")
ax2.set_ylabel("Inferred Elapsed Response Time (Hours)")
ax2.set_title("B. Inferred Redox Signaling Kinetics", weight="bold")
ax2.set_ylim(0, 7.0)
for b in bars:
    y = b.get_height()
    ax2.text(b.get_x() + b.get_width()/2, y + 0.15, f"{y:.2f} h", ha="center", weight="bold")
ax2.legend(loc="upper left")
ax2.grid(axis="y", linestyle="--", alpha=0.5)

fig.suptitle("Figure 3: Arabidopsis ROS Decoder Identifies Antioxidant Shield Failure in Anthocyaninless Fast Plants", fontsize=12, weight="bold", y=0.98)
plt.tight_layout()
fig.savefig(FIG_DIR / "fig3_redox_cvae_decoding.png")
fig.savefig(FIG_DIR / "fig3_redox_cvae_decoding.pdf")
plt.close(fig)

# -------------------------------------------------------------------------
# Figure 4: Radiation Resilience Index (RRI) & Housekeeping Module Collapse
# -------------------------------------------------------------------------
print("Generating Figure 4: Radiation Resilience Index (RRI)...")
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))

# Panel A: RRI Boxplot across Genotype and Radiation Dose
sns.boxplot(data=rri_scores, x="genotype", y="RRI", hue="condition",
            palette={"0_cGy": "#8AB4F8", "40_cGy": "#EA4335"}, ax=ax1, width=0.5)
sns.stripplot(data=rri_scores, x="genotype", y="RRI", hue="condition",
              dodge=True, color="black", alpha=0.6, ax=ax1, size=6)
ax1.set_title("A. Composite Radiation Resilience Index (RRI)", weight="bold")
ax1.set_xlabel("Genotype")
ax1.set_ylabel("Resilience Index (RRI)")
ax1.set_ylim(0.4, 0.9)
ax1.grid(axis="y", linestyle="--", alpha=0.5)

# Fix duplicate legend from stripplot
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles[:2], labels[:2], title="Dose")

# Panel B: Housekeeping Module Preservation Breakdown
sns.barplot(data=rri_scores, x="genotype", y="RRI_module", hue="condition",
            palette={"0_cGy": "#8AB4F8", "40_cGy": "#EA4335"}, ax=ax2, edgecolor="black", capsize=0.1)
ax2.set_title("B. Sustained Housekeeping Module Preservation", weight="bold")
ax2.set_xlabel("Genotype")
ax2.set_ylabel("Preservation Score (r)")
ax2.set_ylim(0, 0.7)
ax2.grid(axis="y", linestyle="--", alpha=0.5)
ax2.legend(title="Dose")

fig.suptitle("Figure 4: 40 cGy GCR Drives RRI Nadir in Wild-Type and Constitutive Stress in Anthocyaninless Plants", fontsize=12, weight="bold", y=0.98)
plt.tight_layout()
fig.savefig(FIG_DIR / "fig4_rri_degradation.png")
fig.savefig(FIG_DIR / "fig4_rri_degradation.pdf")
plt.close(fig)

# -------------------------------------------------------------------------
# Figure 5: Scent Route Metabolic Trade-Off vs Radiation
# -------------------------------------------------------------------------
print("Generating Figure 5: Floral Scent vs Radiation Trade-off...")
fig, ax = plt.subplots(figsize=(9, 5))

routes = scent_tradeoff["route"].str.replace("_", " ").str.title()
wt_deltas = scent_tradeoff["wt_mean_delta"]
al_deltas = scent_tradeoff["al_mean_delta"]

y = np.arange(len(routes))
h = 0.35

ax.barh(y - h/2, wt_deltas, height=h, label="Wild-Type (40 cGy vs 0 cGy)", color="#1E8E3E", edgecolor="black")
ax.barh(y + h/2, al_deltas, height=h, label="Anthocyaninless (40 cGy vs 0 cGy)", color="#D93025", edgecolor="black")

ax.set_yticks(y)
ax.set_yticklabels(routes)
ax.axvline(0, color="black", linestyle="-", linewidth=1)
ax.set_xlabel("Log2 Fold-Change (40 cGy vs 0 cGy)")
ax.set_title("Figure 5: Metabolic Divergence in Floral Scent Pathways Under GCR Radiation", weight="bold")
ax.legend(loc="lower right")
ax.grid(axis="x", linestyle="--", alpha=0.5)

plt.tight_layout()
fig.savefig(FIG_DIR / "fig5_scent_metabolic_tradeoff.png")
fig.savefig(FIG_DIR / "fig5_scent_metabolic_tradeoff.pdf")
plt.close(fig)

# -------------------------------------------------------------------------
# Figure 6: Spatial Anatomical Vulnerability (ggPlantMap)
# -------------------------------------------------------------------------
print("Generating Figure 6: Spatial Anatomical Vulnerability...")
fig, ax = plt.subplots(figsize=(10, 5))

spat_pivot = redox_spat.pivot(index="cell_type", columns="contrast", values="vulnerability_score")
spat_pivot = spat_pivot.rename(columns={
    "WT_GCR_40cGy": "WT GCR 40cGy",
    "AnthoLess_GCR_40cGy": "AnthoLess GCR 40cGy",
    "Interaction_Dose_x_Genotype": "Interaction",
    "Scent_High_vs_Low": "Scent High vs Low"
})
spat_pivot = spat_pivot.sort_values(by="Interaction", ascending=False)

sns.heatmap(spat_pivot * 100, annot=True, fmt=".1f", cmap="YlOrRd", ax=ax, cbar_kws={"label": "Predicted Vulnerability Share (%)"})
ax.set_title("Figure 6: Spatial Anatomical Susceptibility Profile Across Plant Organs", weight="bold")
ax.set_ylabel("Tissue / Cell Layer (ggPlantMap)")

plt.tight_layout()
fig.savefig(FIG_DIR / "fig6_spatial_vulnerability.png")
fig.savefig(FIG_DIR / "fig6_spatial_vulnerability.pdf")
plt.close(fig)

print("=" * 70)
print("20_generate_figures_cross_species: All 6 figures generated successfully!")
print("=" * 70)
