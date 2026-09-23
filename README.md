# Multimodal Physiological Noise Regression in fNIRS: Methodological Validation Through a Novel Finger-Tapping and Breath-Holding Dataset

[![MATLAB](https://img.shields.io/badge/MATLAB-R2022b%2B-blue.svg)](https://www.mathworks.com/products/matlab.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![University of Padova](https://img.shields.io/badge/University%20of%20Padova-DEI-red.svg)](https://www.dei.unipd.it/)
[![Industry Collaboration](https://img.shields.io/badge/Collaboration-NIRx%20Medical%20Technologies-orange.svg)](https://nirx.net/)

Official open-source MATLAB analysis pipeline developed for the Master's Thesis in **Bioengineering for Neuroscience** at the **University of Padova (Department of Information Engineering - DEI)**, conducted in scientific collaboration with **NIRx Medical Technologies LLC** (Berlin, Germany).

---

## 📋 Table of Contents
- [Executive Summary](#executive-summary)
- [Experimental Paradigm & Multimodal Instrumentation](#experimental-paradigm--multimodal-instrumentation)
- [Two-Tier Computational Framework](#two-tier-computational-framework)
- [Hierarchical Model Architectures (M0–M6)](#hierarchical-model-architectures-m0m6)
- [Data Organization & Prerequisites](#data-organization--prerequisites)
- [Pipeline Workflow](#pipeline-workflow)
- [Script Descriptions & Methodology](#script-descriptions--methodology)
- [Key Scientific Findings](#key-scientific-findings)
- [Requirements & Dependencies](#requirements--dependencies)
- [Step-by-Step Execution Guide](#step-by-step-execution-guide)
- [Data Availability & Ethics](#data-availability--ethics)
- [Scientific References](#scientific-references)
- [Citation & Academic Attribution](#citation--academic-attribution)

---

## 🔬 Executive Summary

Systemic physiological fluctuations (including arterial cardiac pulsatility, respiratory chest excursions, blood pressure Mayer waves, hypercapnic cerebral vasodilation, and sympathetic sudomotor bursts) introduce widespread non-neural confounds into functional Near-Infrared Spectroscopy (fNIRS) recordings. When cognitive or motor tasks coincide with systemic autonomic stress, these confounds severely distort hemodynamic response functions (HRFs) and trigger severe false-positive activation explosions across non-motor cortices.

While superficial **short-separation (SC) channels** (d < 15 mm) capture extracerebral scalp shunting, they remain largely blind to systemic vascular surges that propagate through the deeper cerebral arterial tree. This repository provides an end-to-end, leakage-resistant computational pipeline that:
1. Preprocesses and synchronizes multi-channel fNIRS with multi-lead autonomic physiological biosignals (Lead II ECG, earlobe PPG, chest respiration belt, galvanic skin conductance, and pulse oximetry).
2. Empirically proves 1/f^β colored noise properties (β ≈ 1.7) and establishes that block-by-block linear detrending is mathematically mandatory to avoid spurious cross-correlations (|R| > 0.35).
3. Formulates a nested hierarchy of regularized ridge regression models (**M0–M6**) incorporating multi-lag temporal expansions (τ ∈ {0, 2, 4, 6, 8, 10, 12} s).
4. Eliminates data leakage by calibrating regression weights strictly on compliant hypercapnic Breath-Holding (**BH**) blocks via Leave-One-Block-Out Cross-Validation (**LOBO-CV**), transferring frozen mappings to an independent dual-task condition (**COMB**: simultaneous motor execution + breath holding).
5. Benchmarks decontaminated activations against an independent, SC-corrected cortical ground truth estimated via simultaneous deconvolution GLM in Homer2.
6. Evaluates population generalizability across N = 10 participants using paired non-parametric Wilcoxon signed-rank tests with step-down Holm-Bonferroni multiplicity adjustments and the one-standard-error (**1-SE**) parsimony selection rule.

---

## 🧪 Experimental Paradigm & Multimodal Instrumentation

### Participant Cohort
The experimental protocol was conducted on N = 10 healthy adult volunteers (8 males, 2 females; age: mean ± SD = 30.3 ± 5.0 years, range: 25–40 years). All participants were right-handed, ensuring consistent contralateral primary motor cortex activation (C3 region) during unilateral motor execution. All participants provided informed voluntary consent adhering to the Declaration of Helsinki.

### Hardware & Sensor Suite
* **fNIRS Optical Imaging**: Continuous-wave **NIRSport2** system (NIRx Medical Technologies LLC, Berlin, Germany).
  * Dual-wavelength LED illumination at 760 nm and 850 nm.
  * Acquisition rate: 12.60 Hz configured via Smart Spatial Multiplexing with Enhanced Frequency Encoding (EFE) in NIRx Aurora.
  * Probe Layout: **54 optical channels** (16 sources, 15 detectors) mapped onto the international 10–20 system over the motor cortex:
    * **46 long cortical channels** (source–detector separation d ≈ 30 mm).
    * **8 dedicated short-separation channels** (d ≈ 8 mm) co-located with detectors to isolate superficial scalp hemodynamics.
* **Auxiliary Autonomic Sensor**: Wearable **NIRx WINGS2** module streaming 4 channels wirelessly at 500 Hz (24-bit resolution):
  * **Lead II ECG**: Active bipolar ExG with 3 Ag/AgCl disposable electrodes.
  * **Photoplethysmography (PPG)**: Optical earlobe clip sensor measuring peripheral blood volume pulsations.
  * **Respiration Belt**: Piezoelectric/inductive thoracic stretch belt recording chest wall excursions.
  * **Electrodermal Activity (EDA/GSR)**: Constant 0.5 V DC excitation across two palmar finger straps on the left hand.
  * **Pulse Oximetry (SpO2)**: Capillary blood oxygen saturation (70–100%).
* **Synchronization**: Sub-millisecond hardware-level event markers and cross-device alignment established via the **Lab Streaming Layer (LSL)** centralized in NIRx Aurora. Task presentation controlled via **PsychoPy**.

### Experimental Protocol Structure
The randomized block design comprised a total recording duration of 25 minutes (1500 s):
1. **Initial Baseline (30 s)**: Spontaneous resting-state baseline.
2. **Training Phase (4 min)**: 6 practice trials (2 per condition) to verify participant comprehension.
3. **Instruction Pause**: Self-paced rest; participants reviewed on-screen cues instructing end-expiratory apnea initiation.
4. **Main Experiment (20 min, 20 randomized blocks)**:
   * **Finger Tapping (FT, 10 trials)**: Self-paced sequential finger-to-thumb tapping with the right hand (21 s block). Selective contralateral neurovascular coupling in left motor cortex.
   * **Breath Holding (BH, 5 trials)**: End-expiratory voluntary apnea (21 s block). Acute arterial hypercapnia and systemic cerebral vasodilation without motor task execution.
   * **Combined Challenge (COMB, 5 trials)**: Simultaneous right-hand finger tapping + voluntary breath holding (21 s block). Severe stress-test where focal motor activation is directly masked by global hypercapnic vasodilation.
5. **Final Baseline (30 s)**: Terminal resting recovery.
* **Trial Architecture**: 30–42 s pseudorandomized jittered rest → 3 s pre-stimulus visual cue → 21 s active task execution initiated by synchronized LSL trigger (t = 0).

---

## ⚙️ Two-Tier Computational Framework

Signal processing is divided between two complementary environments:

```
[Raw Optical SNIRF & Auxiliary WINGS Files]
                    │
 ┌──────────────────┴──────────────────────────────────────┐
 │  STAGE 1: Optical Preprocessing in NIRx Satori           │
 │  ├── Channel pruning via Coefficient of Variation (CV ≤ 10%)
 │  ├── Optical density transformation: ΔOD(t, λ) = -ln(I/I0)
 │  ├── Monotonic cubic interpolation for transient spike removal
 │  ├── Temporal Derivative Distribution Repair (TDDR) for motion shifts
 │  ├── Zero-phase 2nd-order Butterworth bandpass filter (0.005–0.50 Hz)
 │  ├── Modified Beer-Lambert Law (MBLL) inversion to Δ[HbO] and Δ[HbR]
 │  └── Export preprocessed file: <subject>_Satori2.snirf
 └──────────────────┬──────────────────────────────────────┘
                    ▼
 ┌─────────────────────────────────────────────────────────┐
 │  STAGE 2: Downstream Physiological Modeling in MATLAB   │
 │  ├── physio_general.m: Conditioning, R-peaks, cvxEDA, compliance audit
 │  ├── script_noise_simulation.m: 1/f^β PSD & Monte Carlo detrending test
 │  ├── script_task_correlation.m: Dynamic 8x8 task-state correlation matrices
 │  ├── integrated_spatial_agreement_validation.m: Homer2 GLM & nested ridge LOBO-CV
 │  └── group_analysis_thesis.m: Second-level Wilcoxon tests & 1-SE selection
 └─────────────────────────────────────────────────────────┘
```

---

## 📊 Hierarchical Model Architectures (M0–M6)

To determine whether multimodal sensors provide incremental predictive value beyond superficial optical channels, the pipeline evaluates seven nested model tiers:

| Tier | Architecture | Predictors (P) | Physiological Predictor Components |
|:---:|:---|:---:|:---|
| **M0** | **Raw Baseline** | 0 | Uncorrected concentration data (no regressors) |
| **M1** | **Short Channels (SC)** | 1 | First spatial principal component of SC channels (SC-PCA1) at zero lag (τ = 0) |
| **M2** | **SC + Respiration** | 15 | M1 + Conditioned chest displacement (Resp_raw) & dynamic breathing vigor (Resp_activity) across 7 lags |
| **M3** | **SC + Resp + ECG-HR** | 22 | M2 + Electrophysiological instantaneous heart rate (HR_ECG) across 7 lags |
| **M4** | **Full Cardiorespiratory** | 29 | M3 + Arterial oxygen saturation (SpO2) across 7 lags |
| **M5** | **M4 + Phasic EDA** | 36 | M4 + Sympathetic sudomotor phasic driver (EDA_phasic, SCR from cvxEDA) across 7 lags |
| **M6** | **Full Multimodal** | 43 | M5 + Autonomic background tonic drift (EDA_tonic, SCL from cvxEDA) across 7 lags |

### Multi-Lag Temporal Embedding
Because vascular and autonomic responses propagate with physiological latencies (pulse transit time ≈ 150–250 ms, sympathetic sudomotor delay ≈ 1–3 s, hypercapnic vasodilation ≈ 4–10 s), all peripheral signals in M2–M6 are expanded across seven discrete temporal lags:

$$\tau_j \in \{0, 2, 4, 6, 8, 10, 12\}\text{ seconds}$$

### Regularized Ridge Regression
To resolve ill-conditioning and collinearity among time-lagged physiological regressors, channel-wise parameters are estimated via L2-regularized ridge regression:

$$\hat{\boldsymbol{\beta}}_{i,m}(\lambda) = \left( \mathbf{X}_m^T \mathbf{X}_m + \lambda \mathbf{I}_{P_m} \right)^{-1} \mathbf{X}_m^T \mathbf{y}_i$$

evaluated across an 8-point logarithmic candidate grid:

$$\Lambda = \{0, 10^{-3}, 10^{-2}, 10^{-1}, 1, 10, 100, 1000\}$$

---

## 📁 Data Organization & Prerequisites

Each subject directory within the `DATASET/` folder must contain three raw acquisition files:

```
DATASET/
└── <subject_id>/
    ├── <subject_id>.snirf              # Raw fNIRS file (contains 3D probe geometry: sourcePos3D, detectorPos3D, triggers)
    ├── <subject_id>.wings              # Raw auxiliary physiological recording from NIRx WINGS2
    ├── <subject_id>_Satori2.snirf      # Preprocessed fNIRS file exported from Brain Innovation Satori
    │                                   # (contains motion-corrected, bandpassed, MBLL-inverted Δ[HbO]/Δ[HbR] in μM)
    │
    │   ── Generated sequentially during pipeline execution ──
    ├── <subject_id>_ProcessedPhysio.csv # Output of physio_general.m
    ├── <subject_id>_ValidTrials.mat     # Output of physio_general.m
    ├── <subject_id>_TaskCorrelation.mat # Output of script_task_correlation.m
    ├── <subject_id>_IntegratedSpatial_ModelSummary.csv  # Output of integrated_spatial_agreement_validation.m
    └── <subject_id>_IntegratedSpatial_ActivationAudit.csv # Output of integrated_spatial_agreement_validation.m
```

> **Why both `.snirf` files are required:**  
> `<subject>.snirf` contains the uncalibrated raw optical intensities alongside the physical 3D coordinates of optodes registered to the 10–20 system, required by the scripts to dynamically identify short channels (d < 15 mm) via Euclidean distance norm. `<subject>_Satori2.snirf` contains the preprocessed hemodynamic concentration changes (Δ[HbO] and Δ[HbR] in μM) after TDDR motion correction, bandpass filtering, and MBLL conversion.

---

## 🔄 Pipeline Workflow

```
[Raw Acquisition: <subject>.snirf & <subject>.wings]
           │
           ▼
  1. physio_general.m
     ├── Lead II ECG R-peak detection & PPG peak tracking
     ├── cvxEDA decomposition (Phasic SCR / Tonic SCL)
     ├── Respiration Savitzky-Golay filtering & compliance audit (1.8 * Anorm)
     └── Produces: *_ProcessedPhysio.csv & *_ValidTrials.mat
           │
           │  [Preprocessed fNIRS: <subject>_Satori2.snirf] 
           │  (Satori Stage 1: MBLL Δ[HbO]/Δ[HbR], TDDR, Bandpass 0.005–0.50 Hz)
           │            │
           ├────────────┼─────────────────────────────────┐
           ▼            ▼                                 ▼
  2. script_noise_simulation.m              3. script_task_correlation.m
     ├── 1/f^β PSD characterization (β=1.7)      ├── Dynamic SC PCA extraction (HbO/HbR)
     └── Monte Carlo detrending validation       ├── State-dependent 8x8 correlation matrices
                                                 └── Produces: *_TaskCorrelation.mat
                                                               │
           ┌───────────────────────────────────────────────────┘
           ▼
  4. integrated_spatial_agreement_validation.m
     ├── Empirical FT reference (Simultaneous Homer GLM with SC PCA)
     ├── Nested Ridge LOBO-CV on BH blocks (Models M0-M6)
     ├── Transfer frozen models to COMB challenge
     ├── FDR-corrected activation maps & C3 motor HRF RMSE
     └── Produces: *_IntegratedSpatial_ModelSummary.csv
           │
           ▼
  5. group_analysis_thesis.m
     ├── Paired Wilcoxon signed-rank tests (Holm-Bonferroni corrected)
     ├── 1-SE parsimony model selection rule
     └── Publication-ready 4-panel summary figures & performance tables
```

---

## 📂 Script Descriptions & Methodology

### 1. `physio_general.m`
* **Purpose**: Comprehensive preprocessing of auxiliary physiological signals recorded via NIRx WINGS2 synchronized with fNIRS.
* **Input Files**: `<subject>.wings` and `<subject>.snirf`.
* **Methodology**:
  * **Dynamic Sampling Rate**: Empirically computed from the timestamp gradient (fs ≈ 500 Hz).
  * **Lead II ECG**: Zero-phase 3rd-order Butterworth bandpass (0.5–40.0 Hz), refractory lockout Δt_min = 0.35 s (≤ 171 BPM), adaptive dual-threshold R-peak detection, outlier rejection (40–140 BPM), 5th-order median filter, and PCHIP interpolation to continuous HR_ECG(t).
  * **PPG**: Zero-phase 3rd-order Butterworth bandpass (0.5–5.0 Hz), sliding window (15 s) adaptive peak detection with 0.45 s lockout (≤ 133 BPM), prominence ≥ 0.5 σ_local, height ≥ 0.25 A_max, median filtered and interpolated to HR_PPG(t).
  * **Respiration**: 3rd-order Savitzky-Golay FIR smoothing filter (`sgolayfilt`, frame ≈ 1.5 s) preserving excursion extrema.
  * **Pre-Apnea Compliance Audit**: Evaluates the [-2 s, +2 s] window around task onset. If A_peak > 1.8 · A_norm (95th percentile of resting breathing), the trial is rejected for pre-apnea hyperventilation, setting `tInc_physio = 0` to prevent hyperoxia artifacts.
  * **Respiratory Activity Metric**: Extracts dynamic breathing vigor via sliding-window standard deviation (4 s window), dropping sharply to zero during apnea.
  * **Electrodermal Activity (`cvxEDA`)**: Resampled to 4 Hz, z-scored, and decomposed via convex quadratic programming into sparse sudomotor nerve impulses, phasic driver (τ1 = 0.7 s, τ0 = 2.0 s, α = 0.0008), and smooth cubic B-spline tonic baseline (Δ = 10 s, γ = 0.01).
  * **Pulse Oximetry (SpO2)**: Clamped to 100%, dropouts < 70% imputed via nearest-neighbor interpolation, followed by a 5 s median filter.
* **Outputs**: `<subject>_ProcessedPhysio.csv`, `<subject>_ValidTrials.mat`, Satori `.sdm` design matrix, and multi-panel diagnostic figures.

### 2. `script_noise_simulation.m`
* **Purpose**: Theoretical and empirical validation of 1/f^β low-frequency colored noise in fNIRS and proof of linear detrending necessity.
* **Input Files**: `<subject>.snirf` (for 3D probe geometry), `<subject>_Satori2.snirf`, `<subject>_ProcessedPhysio.csv`.
* **Methodology**:
  * **Periodogram PSD**: Estimates log-log power spectral density of short-channel PC1 (SC-PCA1 for Δ[HbO]) and heart rate (HR_ECG), confirming close adherence to 1/f^1.7 power-law decay (f < 0.10 Hz) with distinct respiratory (≈ 0.28 Hz) and arterial cardiac peaks (≈ 1.24 Hz).
  * **Monte Carlo Simulation (N = 1000 iterations)**: Synthesizes independent colored noise pairs (β = 1.7) across the experimental block timeline. Proves that raw concatenated segments induce severe spurious correlations (|R| ≈ 0.26–0.35).
  * **Detrending Proof**: Demonstrates that block-by-block linear detrending collapses the empirical null distribution to zero (|R| ≤ 0.08–0.09 for FT; |R| ≤ 0.12 for BH/COMB), establishing the empirical significance threshold of |R| ≥ 0.12*.

### 3. `script_task_correlation.m`
* **Purpose**: Maps state-dependent physiological coupling between peripheral autonomic signals and superficial scalp hemodynamics across tasks (extending von Lühmann et al., 2020).
* **Input Files**: `<subject>.snirf`, `<subject>_Satori2.snirf`, `<subject>_ProcessedPhysio.csv`, `<subject>_ValidTrials.mat`.
* **Methodology**:
  * **Dynamic Short-Channel PCA**: Computes Euclidean distance d from 3D coordinates, selects channels with d < 15 mm, rejects inconsistent channels (R < 0.10), and extracts PC1 for Δ[HbO] and Δ[HbR].
  * **Segmented 8x8 Correlation Matrices**: Evaluates cross-modal correlation across four states: Resting Baseline (30 s), Finger Tapping (10 × 21 s), Breath Holding (5 × 21 s), and Combined (5 × 21 s).
  * **Detrending Comparison**: Visualizes raw vs. block-by-block detrended heatmaps, flagging cells exceeding the empirical threshold (|R| ≥ 0.12*).
* **Outputs**: `<subject>_TaskCorrelation.mat` and multimodal time-series overlay figures.

### 4. `integrated_spatial_agreement_validation.m`
* **Purpose**: Subject-level integrated validation pipeline with leakage-resistant architecture.
* **Input Files**: `<subject>.snirf`, `<subject>_Satori2.snirf`, `<subject>_ProcessedPhysio.csv`, `<subject>_ValidTrials.mat`, `<subject>_TaskCorrelation.mat`.
* **Methodology**:
  * **Empirical FT Reference**: Simultaneous deconvolution GLM (FT + BH + COMB) in Homer2 (`hmrDeconvHRF_DriftSS`) using Gaussian basis functions (σ = 0.5 s, Δτ = 0.5 s), 3rd-order polynomial drift, and SC-PCA1 as nuisance regressor. Yields benchmark HRF and FDR-corrected (q < 0.05) binary activation ground truth.
  * **Nested Ridge LOBO-CV on BH**: Evaluates candidate models M0–M6 across compliant BH blocks. Fits shrinkage parameter λ* on N_BH - 1 training blocks and computes nRMSE over post-stimulus window [0, 40 s] across all 46 long channels on the held-out test block.
  * **Model Selection**: Evaluates minimum CV error (m_min) and 1-SE parsimony rule (m_1SE).
  * **Independent COMB Transfer**: Freezes weights from pure hypercapnia and applies them forward to decontaminate COMB: Y_clean,COMB = Y_COMB - X_COMB * B_hat.
  * **Spatial & Temporal Validation**:
    * **Temporal C3 HRF RMSE**: Evaluated over [0, 40 s] in the 12 long channels closest to D4.
    * **Spatial Pearson Correlation r**: Continuous correlation with FT reference map across all 46 channels, with 95% BCa bootstrap confidence intervals (N_boot = 2000).
    * **FDR Activation Maps (q < 0.05)**: Evaluates F1-score and Non-Reference Activation Rate: FP / (FP + TN).
    * **Safety Control**: Neurovascular washout distortion test on unperturbed FT motor data.
* **Outputs**: `<subject>_IntegratedSpatial_ModelSummary.csv`, `<subject>_IntegratedSpatial_ActivationAudit.csv`, C3 temporal HRF plots, and 4-panel scalp topographic maps.

### 5. `group_analysis_thesis.m`
* **Purpose**: Second-level statistical inference treating the **subject** as the independent statistical unit (N = 10).
* **Input**: Directory containing all subject `*_IntegratedSpatial_ModelSummary.csv` tables.
* **Methodology**:
  * **Descriptive Statistics**: Computes Median and Interquartile Range (IQR) for primary (BH nRMSE) and secondary outcomes (C3 HRF RMSE, Spatial r).
  * **Hypothesis Testing**: Performs planned paired two-tailed Wilcoxon signed-rank tests with step-down Holm-Bonferroni multiplicity adjustments (α = 0.05):
    1. M0 vs M1 (SC efficacy)
    2. M1 vs M2 (+ Respiration)
    3. M1 vs M4 (Full Cardiorespiratory)
    4. C3 HRF RMSE: M1 vs M4 (Local Motor Recovery)
    5. Spatial r: M1 vs M4 (Topographical Fidelity)
    6. C3 HRF RMSE: M1 vs Adaptive (1-SE)
    7. Spatial r: M1 vs Adaptive (1-SE)
  * **Parsimony Audit**: Computes frequency distribution of 1-SE selected model architectures across the cohort.
* **Outputs**: `Table1_Group_Model_Performance.csv`, `Table2_Group_Hypothesis_Tests.csv`, `Table3_Model_Selection_Counts.csv`, and publication-quality 4-panel figure.

---

## 🏆 Key Scientific Findings

### 1. Single-Subject Level: Outlier Rescue & False-Positive Elimination
In vascularly reactive participants (exemplified by subject `2026-07-02_003` in Chapter 4):
* **Superficial SC Regression Alone (M1) Triggers Severe Spatial Contamination**: While M1 reduces local scalp variance, it is blind to hypercapnic vasodilation in the cerebral compartment. The GLM misattributes global arterial surges to the motor task regressor, producing a catastrophic **56.4% non-reference false-positive rate** (28/46 channels falsely active across bilateral frontal and parietal cortices).
* **Multimodal Cardiorespiratory Regression (M4) Completely Restores Specificity**: Incorporating respiration, heart rate, and pulse oximetry suppresses all artifactual bilateral channels down to **0.0% false positives**, isolates exactly the 2 true contralateral motor channels surrounding D4/C3, cuts temporal C3 HRF error by **52.4%** (RMSE = 0.756 μM vs. 1.589 μM raw), and achieves optimal spatial correlation with ground truth (r = 0.724, 95% bootstrap CI: [+0.167, +0.788]).
* **Overfitting Beyond M4**: Expanding predictors to electrodermal activity (M5–M6) reintroduces synthetic ringing and worsens cross-validation error due to regressor collinearity.

### 2. Group-Level Level: The Parsimony Rule (N = 10)
Across the complete cohort of N = 10 participants:
* **M1 is Highly Significant Over Raw Baseline**: Transition from M0 to M1 yields a robust, statistically significant reduction in prediction error (p_Holm = 0.0137, median difference = -0.0959).
* **Incremental Multimodal Regressors Show No Population-Wide Superiority**: Adding respiration (M2) or cardiorespiratory features (M4) yields non-significant paired differences at the group level (p_Holm = 1.0000).
* **The 1-SE Parsimony Verdict**: Under the 1-SE parsimony rule on independent breath-holding blocks, **M1 was selected in 50% of subjects (5/10)**, M0 in 2 subjects, M2 in 2 subjects, and M5 in 1 subject (M3, M4, and M6 were never selected).
* **Scientific Conclusion**: In typical participants, superficial short channels capture the dominant shared systemic confound. Therefore, **fixed short-channel regression (M1) represents the most robust, parsimonious group-level baseline**, while **multimodal cardiorespiratory regression (M4) serves as an essential, targeted rescue mechanism for physiologically vulnerable or vascularly hyper-reactive participants**.

---

## 💻 Requirements & Dependencies

### MATLAB Environment
* **MATLAB**: R2022b or later recommended.
* **MathWorks Toolboxes**:
  * Signal Processing Toolbox (`findpeaks`, `butter`, `filtfilt`, `periodogram`, `sgolayfilt`, `resample`, `medfilt1`)
  * Statistics and Machine Learning Toolbox (`pca`, `zscore`, `corr`, `signrank`, `bootstrp`)

### External Open-Source Toolboxes
1. **[Homer2](https://homer-fnirs.org/)**: Required for `hmrDeconvHRF_DriftSS.m` (joint deconvolution GLM with short-channel PCA).  
   *Note: Built-in shadowing protection ensures MATLAB's native `findpeaks` is prioritized over Homer2's legacy version.*
2. **[cvxEDA](https://github.com/lciti/cvxEDA)**: Required for convex optimization decomposition of electrodermal activity (Greco et al., 2016).

---

## 🚀 Step-by-Step Execution Guide

### 1. Installation
Clone this repository into your local MATLAB working directory:
```bash
git clone https://github.com/<your-username>/fnirs-systemic-physiology-regression.git
cd fnirs-systemic-physiology-regression
```

### 2. Add Toolboxes to Search Path
```matlab
addpath(genpath('/path/to/homer2'));
addpath(genpath('/path/to/cvxEDA'));
```

### 3. Run Pipeline Sequentially
```matlab
% Step 1: Preprocess physiological recordings and validate trial compliance
physio_general

% Step 2: Characterize 1/f noise spectrum and run Monte Carlo detrending test
script_noise_simulation

% Step 3: Compute state-dependent multimodal correlation matrices
script_task_correlation

% Step 4: Run single-subject integrated spatial validation (Homer GLM + Ridge LOBO-CV)
results = integrated_spatial_agreement_validation();

% Step 5: Execute group-level inferential analysis across all subjects
group_results = group_analysis_thesis();
```

---

## 🔒 Data Availability & Ethics

Due to European General Data Protection Regulation (GDPR) and ethical commitments regarding human participant neuroimaging, raw fNIRS (`.snirf`), preprocessed hemodynamic files (`*_Satori2.snirf`), and auxiliary physiological recordings (`.wings`) are not publicly archived in this repository. De-identified data matrices and summary statistics are available from the author upon reasonable academic request for non-commercial scientific research.

---

## 📚 Scientific References

1. **von Lühmann, A., et al.** (2020). *Improved physiological noise regression in fNIRS: A multimodal approach.* NeuroImage, 213, 116726.
2. **Tak, S., & Ye, J. C.** (2014). *Statistical analysis of fNIRS data: A comprehensive review.* NeuroImage, 85, 72–91.
3. **Greco, A., et al.** (2016). *cvxEDA: A Convex Optimization Approach to Electrodermal Activity Processing.* IEEE TBME, 63(4), 797–804.
4. **Fishburn, F. A., et al.** (2019). *Temporal Derivative Distribution Repair (TDDR): A motion correction method for fNIRS.* NeuroImage, 184, 171–179.
5. **Guglielmini, M., et al.** (2024). *Hemodynamics and systemic physiology during voluntary breath-holding challenges.* Physiological Measurement, 45, 035002.
6. **Yule, G. U.** (1926). *Why do we sometimes get nonsense-correlations between Time-Series?* Journal of the Royal Statistical Society, 89(1), 1–63.
7. **Zarahn, E., et al.** (1997). *Empirical evaluation of 1/f noise in fMRI and its statistical implications.* NeuroImage, 5(3), 179–197.
8. **Schafer, R. W.** (2011). *What Is a Savitzky-Golay Filter?* IEEE Signal Processing Magazine, 28(4), 111–117.

---

## 👩‍🔬 Citation & Academic Attribution

If you utilize this pipeline, scripts, or methodology in your research, please cite the associated Master's Thesis:

```bibtex
@mastersthesis{monego2026multimodal,
  author       = {Silvia Monego},
  title        = {{Multimodal Physiological Noise Regression in fNIRS: Methodological Validation Through a Novel Finger-Tapping and Breath-Holding Dataset}},
  school       = {University of Padova},
  department   = {Department of Information Engineering (DEI)},
  year         = {2026},
  type         = {Master's Thesis in Bioengineering for Neuroscience},
  address      = {Padova, Italy},
  note         = {In collaboration with NIRx Medical Technologies LLC, Berlin, Germany}
}
```

**Author**: Silvia Monego  
**Degree**: Master of Science in Bioengineering for Neuroscience  
**Institution**: Department of Information Engineering (DEI), University of Padova, Italy  
**Academic Year**: 2025/2026  
