# PGA-CFAR: Polarimetric Sea-State Aware Ship Detection

## 1. Overview

This repository implements a Polarimetric Sea-State Aware Constant False Alarm Rate (PGA-CFAR) framework for ship detection in ocean Synthetic Aperture Radar (SAR) imagery.

The proposed method addresses the limitations of conventional CFAR detectors in heterogeneous ocean clutter by incorporating polarimetric information and sea-state awareness. The system integrates feature extraction, sea-state modelling, and adaptive detection to improve robustness and reduce false alarms.

---

## 2. Objectives

* To model ocean clutter variability using polarimetric features
* To design an adaptive CFAR detector guided by sea-state conditions
* To reduce false alarms in heterogeneous ocean environments
* To improve detection performance near clutter transitions

---

## 3. Key Components of the Method

### 3.1 Polarimetric Feature Extraction

From the coherency matrix, the following features are computed:

* Entropy (H)
* Mean Alpha (α)
* Anisotropy (A)
* Span power
* Co-polarization Phase Difference (CPD)
* Bragg ratio

These features capture the physical scattering behavior of ocean surfaces and ship targets.

---

### 3.2 Polarimetric Sea-State Index (PSSI)

A composite index is defined as:

PSSI = w1·H + w2·Φ + w3·ΔB

where:

* H represents scattering randomness
* Φ represents CPD deviation
* ΔB represents Bragg ratio deviation

The PSSI provides a spatial map of sea-state variability.

---

### 3.3 Adaptive CFAR Detection

The detection stage consists of:

* Adaptive threshold scaling using PSSI
* Elliptical background window aligned with wave direction
* Polarimetric pre-screening using H and α

This allows the detector to adapt dynamically to local clutter conditions.

---

## 4. Repository Structure and File Description

### 4.1 Main Execution Script

* **MAIN_PGA_CFAR.m**

  * Entry point of the project
  * Executes the full pipeline sequentially
  * Generates feature maps, detection outputs, and evaluation results

---

### 4.2 Module Descriptions

* **module1_polsar_reader.m**

  * Loads quad-polarimetric SAR data
  * Supports multiple formats (RADARSAT-2, GF-3, ALOS-2)
  * Generates synthetic scene if real data is unavailable
  * Outputs: SHH, SHV, SVV, incidence angle map

* **module2_build_coherency.m**

  * Constructs the 3×3 coherency matrix using Pauli decomposition
  * Applies spatial averaging for stability

* **module2_lee_filter_T.m**

  * Implements refined Lee filter
  * Reduces speckle noise while preserving polarimetric structure

* **module3_feature_engine.m**

  * Computes all polarimetric features:

    * Entropy (H)
    * Alpha (α)
    * Anisotropy (A)
    * CPD
    * Bragg ratio
    * Span

* **module4_pssi_mapper.m**

  * Computes PSSI from extracted features
  * Applies Gaussian smoothing
  * Computes gradient magnitude and direction

* **module5_pga_cfar.m**

  * Implements PGA-CFAR detection algorithm
  * Includes:

    * Adaptive threshold scaling
    * Elliptical kernel construction
    * Polarimetric pre-screen
  * Also contains baseline implementations:

    * CA-CFAR
    * OS-CFAR
    * H-CFAR

* **module6_post_processor.m**

  * Performs post-processing:

    * Morphological filtering
    * Connected component labeling
    * Centroid extraction
  * Generates final ship detections

---

## 5. Execution Instructions

### 5.1 Requirements

* MATLAB (R2020 or later recommended)
* Statistics Toolbox (optional; only required for gamma random generation)

---

### 5.2 Steps to Run

1. Open MATLAB
2. Set the working directory to the project folder
3. Run the main script:

```matlab
MAIN_PGA_CFAR
```

4. The pipeline will execute all modules sequentially

---

### 5.3 Output Location

Results are saved in the `results/` directory and include:

* Detection maps
* Feature maps
* ROC curves
* CSV summaries

---

## 6. Output Files and Their Meaning

* **detection_results.mat**

  * Contains detected ship centroids and cluster areas
  * Represents final output of the pipeline

* **detection_summary.csv**

  * Tabular summary of detected ships
  * Includes coordinates and sizes

* **polarimetric_features.mat**

  * Stores all intermediate feature maps
  * Useful for analysis and debugging

* **pssi_statistics.csv**

  * Provides statistical distribution of PSSI
  * Separates ocean and ship pixel behaviour

* **roc_data.mat**

  * Contains detection performance metrics across thresholds
  * Used to generate ROC curves

---

## 7. Interpretation of Results

* Ship targets exhibit:

  * High alpha values (double-bounce scattering)
  * Low entropy compared to rough ocean
  * Distinct clustering in feature space

* PSSI map:

  * Low values correspond to calm ocean
  * High values correspond to rough sea
  * Ship pixels remain distinguishable due to entropy suppression

* Detection performance:

  * Reduced false alarms compared to conventional CFAR
  * Improved detection near clutter boundaries
  * Strong ROC performance indicating robustness

---

## 8. Baseline Comparison

The repository includes comparison with:

* CA-CFAR (standard method)
* OS-CFAR (order-statistic based)
* H-CFAR (entropy-adaptive)

The proposed PGA-CFAR demonstrates:

* Higher detection accuracy
* Lower false alarm rates
* Better adaptability to heterogeneous clutter

---

## 9. Authors

* Amritha S
* Yugeshwaran P

---

## 10. Repository Link

https://github.com/Amritha902/PGA-CFAR

---

## 11. Notes

* Synthetic dataset is used for controlled validation
* The framework is modular and extensible
* Suitable for research and academic applications in SAR processing

---
