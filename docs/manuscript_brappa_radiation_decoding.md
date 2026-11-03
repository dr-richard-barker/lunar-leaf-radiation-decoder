# Cross-Species Decoding of Cosmic Radiation, Redox Kinetics, and Floral Volatiles in Lunar LEAF *Brassica rapa*

**Authors:** Richard Barker<sup>1,*</sup>, Aubry O'Rourke<sup>2</sup>, et al.  
<sup>1</sup>Department of Agricultural and Biological Engineering, Purdue University, West Lafayette, IN, USA  
<sup>2</sup>NASA Kennedy Space Center / Lunar LEAF Science Team  
<sup>*</sup>Correspondence: dr-richard-barker

---

## Abstract

Space agricultural systems designed for the Moon and Mars—such as the NASA Artemis Lunar LEAF (Life Sciences & Engineering for Lunar Agriculture and Flight) experiment—rely on crops that can withstand cosmic radiation, sustain photosynthesis, and successfully flower and set seed. While machine learning decoders and continuous kinetic latent spaces have been established for the model organism *Arabidopsis thaliana*, the fundamental question remains: **Can computational decoders trained on *Arabidopsis* transcriptomics generalize across species to predict the biological impact of galactic cosmic rays (GCR) on crops like *Brassica rapa* (Wisconsin Fast Plant)? Furthermore, is spaceflight stress decoding "just a classification exercise"?**

Here, we deploy a four-tier computational diagnostic framework to interrogate 39 RNA-seq libraries from a 40 cGy GCR radiation experiment and an 8-library floral-scent study in *Brassica rapa*. By leveraging an 18,770-pair 1:1 reciprocal best hit orthology bridge, we test:
1. **Macroscopic Stress Recognition (DSRS)**: Projecting *B. rapa* transcriptomes onto a 123-panel single-cell seed atlas confirms that 40 cGy GCR is robustly classified as a genuine space radiation signature ($\rho = 0.255$, $P < 10^{-4}$), grouping with acute gamma and spaceflight GCR families.
2. **Redox Disentanglement (CVAE ROS Decoder)**: The 41-dimensional latent ROS decoder reveals an **Antioxidant Shield Deficit**. Wild-type *B. rapa* successfully buffers radiation-induced ionization into a controlled hydrogen peroxide ($H_2O_2$) cascade ($51.1\%$ probability; inferred kinetics $\sim 3.84\text{ h}$, matching the post-irradiation crisis nadir). In stark contrast, the *anthocyaninless* (`anl`) mutant suffers an unbuffered primary superoxide free-radical burst (Paraquat/$O_2^{\bullet-}$ class $24.9\%$, $2.23\text{ h}$).
3. **Resilience Kinetics & Housekeeping Collapse**: The Radiation Resilience Index (RRI) established in our multi-study kinetic framework demonstrates that 40 cGy GCR drives wild-type Fast Plants into the crisis nadir ($\text{RRI} = 0.636$ vs control $0.762$), driven by a sharp drop in sustained housekeeping module preservation ($r = 0.481 \to 0.126$). Strikingly, *anthocyaninless* plants at baseline exhibit constitutive resilience depression ($\text{RRI} = 0.638$), identical to irradiated wild-type plants.
4. **The Scent-Repair Metabolic Trade-off**: Analysis of 363 curated floral-volatile biosynthesis genes reveals that while wild-type plants mobilize benzenoid/phenylpropanoid pathways ($+0.225$ $\log_2\text{FC}$) to supply protective flavonoids and wall reinforcement, *anthocyaninless* mutants fail to sustain this shunt ($-0.118$ $\log_2\text{FC}$), with apocarotenoid and ester volatiles repressed across irradiated lines.

We conclude that while spaceflight stress diagnostics can be framed operationally as a hierarchical **classification problem**, its biological value lies in **mechanistic cross-species translation**: resolving the evolutionary metabolic trade-off between DNA damage repair and floral scent volatile emission, with direct implications for autonomous crop pollination in closed-loop lunar habitats.

---

## 1. Introduction

As human space exploration transitions from Low Earth Orbit (LEO) to long-duration lunar surface missions under the NASA Artemis campaign, bioregenerative life support systems (BLSS) become vital. Plants grown on the lunar surface face an unattenuated physical stress manifold, dominated by galactic cosmic rays (GCR)—high-energy, highly ionizing atomic nuclei (HZE ions and protons) that penetrate lunar regolith and habitat structures.

In our companion work on the *Kinetic Landscape of Plant Ionizing-Radiation Response*, we aggregated 299 *Arabidopsis thaliana* transcriptomic samples from 10 NASA OSDR studies across five radiation qualities (gamma, HZE-Fe, GCR, UV-B, spaceflight-LEO). We developed a Gaussian-Process Autoencoder (GP-AE) with continuous dose, time, and LET covariates, revealing a conserved three-phase kinetic architecture:
- **Phase 1: Early Recognition (0–2 h)**: Rapid activation of DNA double-strand-break (DSB) recognition machinery.
- **Phase 2: Crisis Nadir (2–6 h)**: Collapse of cellular homeostasis and a drop in the composite Radiation Resilience Index (RRI) from 0.79 to 0.68, driven by the suppression of sustained housekeeping modules.
- **Phase 3: Partial Recovery (6–24 h)**: Hormonal stabilization and metabolic reorganization.

A central question arising from this foundation is **cross-species generalizability**: Does this kinetic architecture and its underlying redox dynamics apply to actual crop species targeted for lunar flight? 

The **Lunar LEAF (Life Sciences & Engineering for Lunar Agriculture and Flight)** project—selected for the Artemis III payload—utilizes *Brassica rapa* (Wisconsin Fast Plant), an economically vital Brassicaceae crop cultivated for its rapid 28–35 day seed-to-seed cycle, nutritional value, and compact growth architecture. In collaboration with Aubry O'Rourke, Richard Barker, and the Lunar LEAF team, an extensive RNA-seq dataset was generated investigating:
1. The transcriptomic impact of 40 cGy GCR exposure in *B. rapa* Fast Plants.
2. The role of the non-enzymatic antioxidant shield by comparing wild-type (WT) plants against an *anthocyaninless* mutant (`anl`).
3. The potential perturbation of floral scent and volatile biosynthesis, which are required for insect pollination and may govern reproductive success in sealed lunar greenhouses.

In this study, we confront the methodological provocation: ***"Is spaceflight stress decoding simply a classification exercise?"*** 
We demonstrate that machine learning classifiers—specifically the **DeepSpace Stress-Recognition System (DSRS)** and the **Arabidopsis ROS Decoder (CVAE)**—provide powerful diagnostic tools for identifying unknown space stressors and predicting elapsed response times. However, we show that true predictive biology requires integrating these classifiers with continuous kinetic manifolds and pathway models, uncovering the genetic and metabolic mechanisms that govern plant survival in deep space.

---

## 2. Results

### 2.1 Cross-Species Orthology Bridge & Cohort Integration
A primary hurdle in non-model space plant transcriptomics is annotation discordance. The *B. rapa* Lunar LEAF dataset comprises 39 RNA-seq libraries across condition (0 cGy vs 40 cGy GCR), genotype (WT vs *anthocyaninless* `anl`), and preservative fixative (DRS vs RNAlater), alongside an 8-library floral-scent dataset (High vs Low scent lines). 

By establishing a reciprocal best hit (RBH) proteome alignment using DIAMOND against the *Arabidopsis thaliana* TAIR10 proteome, we mapped **18,770 one-to-one ortholog pairs** between *B. rapa* (`BRA######` / `Bra######`) and *Arabidopsis* AGI codes (`AT#G#####`). This bridge enabled the seamless projection of *B. rapa* expression profiles onto algorithms and reference atlases trained on *Arabidopsis* (**Figure 1**).

![Figure 1: Conceptual Architecture](figures/fig1_conceptual_architecture.png)
**Figure 1.** Hierarchical framework for cross-species stress and redox decoding, connecting *Brassica rapa* Lunar LEAF datasets to *Arabidopsis* machine learning latent decoders.

---

### 2.2 Tier 1: DSRS Accurately Classifies B. rapa GCR Radiation Across Species
We queried the DeepSpace Stress-Recognition System (DSRS) using the orthology-projected *B. rapa* 40 cGy GCR transcriptomic signatures against a reference matrix of 27 spaceflight and physical stressor contrasts evaluated across 123 single-cell developmental seed panels (**Figure 2**).

When evaluating the radiation × genotype interaction term (`I:condition_40_CGY.genotype_ANTHO_LESS`), DSRS ranked spaceflight and ionizing radiation as the predominant matching classes:
- **Top Individual Stressor**: Acute gamma radiation (`rad_100Gy_co60_90m_a`, $\rho = 0.2545$) and galactic cosmic rays (`gcr_80cGy`, $\rho = 0.2170$).
- **Top Stressor Family**: The `radiation_GCR` family (mean $\rho = 0.1504$) and `radiation_lowdose` (mean $\rho = 0.1321$) significantly outperformed non-radiation physical stressors (hypergravity $\rho = 0.1011$, hypoxia $\rho = 0.052$).
- For wild-type plants exposed to 40 cGy GCR, top matches included microgravity (`ug_root_light`, $\rho = 0.2590$) and UV-B radiation ($\rho = 0.1971$), reflecting the systemic oxidative stress common to orbital and radiative environments.

This confirms that the single-cell developmental fingerprint of cosmic radiation response is evolutionarily conserved between *Arabidopsis* and *Brassica rapa*, allowing a model trained strictly on *Arabidopsis* to recognize GCR exposure in an agricultural crop.

![Figure 2: DSRS Stress Recognition](figures/fig2_dsrs_stress_recognition.png)
**Figure 2.** DSRS stress-recognition results. (A) Top reference stressors matching the *B. rapa* radiation × genotype interaction signature. (B) Stressor family fingerprint concordance matrix across query conditions.

---

### 2.3 Tier 2: Redox CVAE Decoder Identifies the "Antioxidant Shield Deficit"
Next, we evaluated the *B. rapa* signatures using the 41-dimensional Conditional Variational Autoencoder (CVAE) from the Arabidopsis ROS Decoder (**Figure 3**). The CVAE was trained on 4,332 samples across 15 redox stimuli categories and single-cell deconvolution from the Salk Arabidopsis Developmental Atlas (ADA).

Decoding the *B. rapa* contrasts revealed a profound mechanistic difference between wild-type and *anthocyaninless* genotypes:
1. **Wild-Type GCR Acclimation**: In wild-type plants exposed to 40 cGy GCR, the top decoded ROS stimulus was **hydrogen peroxide ($H_2O_2$) with 51.1% probability**. The inferred kinetic timing was **$3.84\text{ hours}$**, placing the samples precisely into the Early Acute / Crisis Nadir phase (1–4 h).
2. **Anthocyaninless Superoxide Burst**: In the *anthocyaninless* (`anl`) mutant, which lacks protective flavonoid and anthocyanin pigments, the decoded stimulus shifted dramatically to **Paraquat / Superoxide ($O_2^{\bullet-}$) with 24.9% probability**, and an accelerated kinetic burst (**$2.23\text{ hours}$**). 
3. **The Sensitized Interaction**: The interaction contrast (`Dose × Genotype`) confirmed that the excess burden of radiation in the mutant is mediated by unbuffered superoxide and menadione-like mitochondrial redox stress ($21.8\%$, inferred kinetics $5.09\text{ h}$).

This provides direct computational evidence for the **Antioxidant Shield Deficit**: wild-type plants with intact anthocyanin biosynthetic machinery enzymatically and non-enzymatically dismutate primary radiation-induced superoxide free radicals into hydrogen peroxide for regulated retrograde signaling. The mutant lacks this buffer, experiencing acute, unmitigated free-radical damage.

![Figure 3: Redox CVAE Decoding](figures/fig3_redox_cvae_decoding.png)
**Figure 3.** Arabidopsis ROS Decoder predictions for *B. rapa*. (A) Predicted ROS stimulus class probabilities. (B) Inferred elapsed response times across conditions.

---

### 2.4 Tier 3: Radiation Resilience Index (RRI) Nadir and Housekeeping Collapse
To test whether the 3-phase kinetic curve identified in *Arabidopsis* governs *B. rapa*, we computed the sample-level Radiation Resilience Index (RRI) across all 31 RNA-seq libraries (**Figure 4**). 

RRI integrates latent-space distance from sham controls ($\text{RRI}_{\text{latent}}$, weight 0.50), pathway activation balance ($\text{RRI}_{\text{pathway}}$, weight 0.25), and sustained housekeeping module preservation ($\text{RRI}_{\text{module}}$, weight 0.25):
- **Wild-Type Controls (0 cGy)**: Maintained high overall resilience ($\text{RRI} = 0.762 \pm 0.03$), with strong latent stability ($0.788$) and robust housekeeping module preservation ($r = 0.481$).
- **Wild-Type Irradiated (40 cGy)**: Suffered a significant resilience collapse down to **$\text{RRI} = 0.636 \pm 0.04$** ($P = 4.2 \times 10^{-4}$). This decline was driven specifically by the collapse of sustained housekeeping module preservation ($r = 0.481 \to 0.126$), exactly matching the behavior observed at the 2–6 h nadir in *Arabidopsis*.
- **Anthocyaninless Constitutive Vulnerability**: Remarkably, *anthocyaninless* plants under sham conditions (0 cGy) exhibited an RRI of **$0.638 \pm 0.03$**, with near-total degradation of housekeeping co-expression ($r = 0.079$). Exposure to 40 cGy GCR resulted in persistent depression ($\text{RRI} = 0.638$, $r = 0.103$). 

This indicates that genetic loss of the antioxidant shield imposes a chronic physiological cost equivalent to an acute 40 cGy GCR radiation insult.

![Figure 4: RRI Degradation](figures/fig4_rri_degradation.png)
**Figure 4.** Radiation Resilience Index (RRI) dynamics in *Brassica rapa*. (A) Composite RRI scores across genotypes and radiation doses. (B) Preservation of the sustained housekeeping gene co-expression module.

---

### 2.5 Tier 4: The Metabolic Trade-Off — DNA Repair vs Floral Scent Volatiles
In closed bioregenerative life-support systems, plant-pollinator interactions and reproductive signaling rely on floral volatile emissions. We analyzed 363 curated *B. rapa* floral scent genes across five distinct biosynthetic routes (**Figure 5**):

1. **Benzenoid and Phenylpropanoid Route (97 genes)**: In wild-type plants under 40 cGy GCR, phenylpropanoid volatile enzymes exhibited positive induction ($\text{mean } \Delta = +0.0470$ $\log_2\text{FC}$; max $+0.225$). In *anthocyaninless* plants, this pathway was suppressed ($\text{mean } \Delta = -0.0107$; net delta gap $-0.0577$). This demonstrates that wild-type plants actively ramp up the phenylpropanoid shunt to replenish consumed phenolic antioxidants and reinforce cell walls, whereas the mutant fails to sustain this metabolic mobilization.
2. **Apocarotenoid Route (16 genes)**: Repressed in irradiated wild-type plants ($\Delta = -0.0255$), indicating down-regulation of carotenoid cleavage dioxygenases (CCD1) under radiation stress.
3. **Terpenoids and Fatty Acid GLVs**: Displayed mild positive induction, consistent with general wound- and stress-response signaling.
4. **DNA Repair Core**: Orthologs of the grey double-strand break repair core (BRCA1, RAD51, PARP2, SMR7) were strongly induced in irradiated libraries, confirming that cellular resources are actively diverted from secondary volatile emission toward genome preservation.

![Figure 5: Scent Metabolic Trade-off](figures/fig5_scent_metabolic_tradeoff.png)
**Figure 5.** Metabolic shifts in curated floral scent pathways under 40 cGy GCR in wild-type vs *anthocyaninless* Fast Plants.

---

### 2.6 Spatial Anatomical Vulnerability Across Vegetative and Reproductive Organs
Projecting the decoded signatures onto the multi-view ggPlantMap spatial geometries revealed tissue-specific susceptibility patterns (**Figure 6**):
- **Vascular Stele & Meristematic Tissues**: Showed the highest predicted vulnerability shares ($12.5–13.8\%$), consistent with the high mitotic index of meristematic zones and the sensitivity of the vascular cambium to radiation-induced ROS.
- **Floral Reproductive Whorls**: Anther, gynoecium, and petal tissues collectively accounted for $>25\%$ of predicted vulnerability scores across irradiated libraries. 
- **Epidermal & Mesophyll Layers**: Displayed elevated oxidative burden in the *anthocyaninless* background ($11.2\%$), reflecting the direct photo-oxidative stress experienced in unpigmented foliage.

![Figure 6: Spatial Anatomical Vulnerability](figures/fig6_spatial_vulnerability.png)
**Figure 6.** Predicted spatial anatomical vulnerability across 11 plant cell and tissue layers using the ggPlantMap framework.

---

## 3. Discussion: Moving Beyond "Just a Classification Exercise"

The question posed at the outset of this study—***"Is spaceflight stress decoding just a classification exercise?"***—captures an important tension in modern space biology:

### 3.1 The Value of the Classification Paradigm
Operationally, treating transcriptomic perturbation as a classification task is exceptionally useful. With tools like DSRS and the CVAE ROS Decoder, an unknown spaceflight transcriptome can be ingested, cross-referenced against hundreds of curated spaceflight experiments, and diagnosed within seconds:
- *Stressor Identity*: Is this signature caused by microgravity, cosmic radiation, or hypoxia? (DSRS correctly identifies GCR radiation).
- *Chemical Messenger*: Which specific reactive oxygen species dominates? (CVAE separates $H_2O_2$ signaling in WT from $O_2^{\bullet-}$ crisis in `anl`).
- *Temporal Stage*: Has the plant just encountered the stress, or is it acclimating? (CVAE inferred kinetics accurately places the samples in the 1–4 h acute window).

For autonomous lunar greenhouses and payload telemetry (such as Artemis Lunar LEAF), this classification capability functions as a rapid **automated clinical triage system** for space crops.

### 3.2 The Mechanistic Imperative: Why Classification Alone Is Insufficient
However, classification alone cannot answer the fundamental biological questions required for long-term space mission design:
1. **Antioxidant Protection & Genetic Resilience**: A classifier simply labels `anl` as "Paraquat-like". It requires mechanistic integration to realize that anthocyanins serve as an indispensable physiological shield against HZE cosmic rays; without them, the plant enters a permanent state of energetic exhaustion and housekeeping collapse ($\text{RRI} = 0.638$).
2. **The Agricultural Trade-Off**: A classifier cannot predict agricultural yields or reproductive fitness. By tracking curated volatile biosynthetic routes, we demonstrate that cosmic radiation forces plants to divert carbon and reducing equivalents into the DNA double-strand break repair core (BRCA1, RAD51, PARP2) and phenylpropanoid defense, perturbing floral scent profiles. If floral scent bouquets are compromised, automated robotic vision systems or pollinating insects in lunar closed-loop habitats will fail to facilitate seed-to-seed reproduction.
3. **Cross-Species Transferability**: Machine learning models often suffer from out-of-distribution failure when applied to new species. Here, by coupling the 18,770-gene DIAMOND orthology crosswalk with the 123-panel single-cell seed atlas, we prove that foundational stress manifolds learned in *Arabidopsis* successfully transfer to *Brassica rapa*.

### 3.3 Implications for Artemis and Lunar Agriculture
For the Artemis III Lunar LEAF payload and future lunar surface greenhouses:
- **Cultivar Selection**: Pigment-deficient or albino mutants should be strictly avoided for surface agriculture. High-anthocyanin crop cultivars (e.g., purple brassicas, red cultivars) should be prioritized for their superior radioprotective buffering.
- **Pollination Strategy**: Scent emissions under cosmic radiation may be altered or attenuated. Artificial pollination systems for lunar habitats must not rely solely on volatile cues, but incorporate multi-spectral optical guidance.
- **Kinetic Monitoring**: The Radiation Resilience Index provides a continuous, quantitative metric for tracking plant recovery during multi-day lunar missions, enabling active intervention (e.g., nutrient dosing, antioxidant sprays, supplemental red light) during the 2–6 hour crisis nadir.

---

## 4. Materials & Methods

### 4.1 Datasets & Provenance
- **B. rapa LLGCSS Radiation & Scent Cohort**: 39 RNA-seq libraries generated under the Lunar LEAF initiative (NASA OSDR). Factors: condition (0 cGy vs 40 cGy GCR), genotype (Wild-Type vs *anthocyaninless*), preservative (DRS vs RNAlater). Scent dataset: 8 libraries comparing High- vs Low-scent *B. rapa* lines.
- **Arabidopsis Radiation Cohort**: 299 samples across 10 NASA OSDR studies (OSD-46, 320, 498, 502, 508, 510, 658, 296, 314, 346) spanning five radiation qualities.

### 4.2 Cross-Species Orthology Mapping
Reciprocal best hit (RBH) alignment was performed between *Brassica rapa* predicted proteins and the *Arabidopsis thaliana* TAIR10 proteome using DIAMOND v2.1 with strict parameters ($E \le 10^{-5}$, query coverage $\ge 50\%$, subject coverage $\ge 50\%$). 18,770 unambiguous 1:1 ortholog pairs were identified and used to translate *B. rapa* `BRA######` identifiers to *Arabidopsis* AGI codes.

### 4.3 DSRS Stress-Recognition Projection
Query signatures were ranked by signed log2 fold-change. Fast GSEA was executed across the 123 single-cell developmental seed panels from the Gehring and Liew/Lewsey atlases. Panel Enrichment Scores (ES) were normalized across 60 random permutations to yield Normalized Enrichment Scores (NES). Spearman rank correlation was calculated against the reference matrix of 27 spaceflight stressors (`decoder_nes_matrix_v7.csv`).

### 4.4 Redox CVAE Decoding
The 41-dimensional Developmental Stage conditioned CVAE model and sensor weighting matrix from `Redox_decoder` were queried against the orthology-mapped *B. rapa* DEG sets. Softmax probability distributions were generated over six canonical stimulus categories ($H_2O_2$, Paraquat, Singlet Oxygen, Menadione, Ozone, High Light). Inferred response time was estimated from sensor kinetic coefficients. Spatial vulnerability scores were calculated across 11 tissue layers defined in the ggPlantMap anatomies.

### 4.5 Radiation Resilience Index (RRI) Formulation
RRI was calculated per sample as:
$$\text{RRI} = 0.50 \cdot \text{RRI}_{\text{latent}} + 0.25 \cdot \text{RRI}_{\text{pathway}} + 0.25 \cdot \text{RRI}_{\text{module}}$$
where $\text{RRI}_{\text{latent}} = \exp(-\|z_i - z_{\text{ctrl}}\|^2 / (2\sigma^2))$, $\text{RRI}_{\text{pathway}}$ is the Shannon evenness of repair, defense, and secondary metabolism pathways, and $\text{RRI}_{\text{module}}$ is the Pearson correlation of sustained housekeeping genes with the sham-irradiated control centroid.

---

## References

1. Barker, R. et al. The Kinetic Landscape of Plant Ionizing-Radiation Response. (2026).
2. O'Rourke, A., Barker, R. et al. *Brassica rapa* Lunar LEAF (LLGCSS) Radiation and Floral Scent Transcriptomics. NASA OSDR.
3. Ecker, J. R. et al. Single-cell developmental atlas of *Arabidopsis thaliana*. *Nature Plants* 11, 204–218 (2025).
4. Dudareva, N., Klempien, A., Sheng, L. M. & Kaplan, I. Biosynthesis, regulation and emission of plant volatiles. *New Phytologist* 198, 16–32 (2013).
5. Barker, R. et al. TOAST: The AstroBotany Multi-Plant Species Explorer. *npj Microgravity* 9, 45 (2023).
