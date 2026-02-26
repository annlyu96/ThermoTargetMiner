
# ThermoTargetMiner

**ThermoTargetMiner** is an online resource for exploring drug-protein interactions in lung cancer research. It leverages the **Proteome Integral Solubility Alteration (PISA)** assay to provide a comprehensive database of targets for 67 compounds.

## Quick Start
Access the live app here: [ThermoTargetMiner Web Tool](https://thermotargetminer.serve.scilifelab.se/app/thermotargetminer)

---
## User Guide

### 1. Sidebar Configuration
- **Dataset Selection**: Select **Cell Line** (A549/NCI-H82) and **Sample Type** (Lysate/In-cell).
- **Drug Selection**: Choose from 67 available drugs.
- **Threshold**: Set the significance level (Default: **4.5** on the `-log10(p)` scale) to identify "pro-targets."

### 2. Data Analysis & Visualization
The right-hand panel displays results in the following order:

#### I. OPLS-DA Model Data
View the summary statistics of the OPLS-DA model for the selected drug to ensure model quality and variance explanation.

#### II. Interactive Scatter Plots (OPLS & Transformed OPLS)
These plots display the protein solubility shifts:
- **Interactive Hover**: Hover over dots to identify specific proteins.
- **Deep-Dive Click**: Click on any protein to see its behavior (solubility levels) across all **four experimental conditions** (A549/H82 x Lysate/In-cell) simultaneously.

#### III. GO Pathway Enrichment
Explore the biological functions of the identified pro-targets. This section highlights enriched Gene Ontology terms to reveal the drug's mechanism of action.

#### IV. Pro-target Count & Overlap
Review the number of significant targets identified and their distribution/overlap across the different experimental setups.

#### V. Pro-target Table
A comprehensive table for browsing and exporting the list of candidate targets.

---
## Scientific Background
The app uses OPLS-DA to "purify" the drug-induced solubility signals from background noise.
- **Lysate Data**: Primarily reveals direct drug-protein binding.
- **In-cell Data**: Reveals the "cellular footprint," including metabolic effects and downstream signaling.
- **Validation**: A protein consistently appearing across multiple datasets (e.g., high k-level) is a high-confidence target.

---
## 📝 Citation
Lyu, H., et al. "ThermoTargetMiner as a proteome integral solubility alteration target database for prospective drugs against lung cancer."
---


