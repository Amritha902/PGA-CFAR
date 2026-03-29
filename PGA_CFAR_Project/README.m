% ==========================================================================
%  PGA-CFAR Ship Detection Framework — README
%  Group 11 | VIT Chennai | ECE B.Tech
%  Course: Radar Signal Processing (Phase 1 → Phase 2)
% ==========================================================================
%
%  HOW TO RUN
%  ----------
%  1. Open MATLAB (R2020b or later recommended).
%     Required toolboxes: Image Processing, Statistics & Machine Learning.
%
%  2. Set MATLAB's Current Folder to the root of this project:
%       >> cd('path/to/PGA_CFAR_Project')
%
%  3. (SYNTHETIC / TESTING MODE — no real data needed)
%     Leave DATA_FORMAT = 'SYNTHETIC' in MAIN_PGA_CFAR.m and run:
%       >> MAIN_PGA_CFAR
%
%  4. (REAL DATA MODE)
%     a. Copy your quad-pol SAR data files into the ./data/ folder.
%     b. Set DATA_FORMAT in MAIN_PGA_CFAR.m:
%          'RADARSAT2' — imagery_HH.bin, imagery_HV.bin,
%                        imagery_VH.bin, imagery_VV.bin, product.xml
%          'ALOS2'     — IMG-HH-*, IMG-HV-*, IMG-VV-* (CEOS format)
%          'GAOFEN3'   — at least 3 GeoTIFF .tif files
%     c. (Optional) Provide an AIS CSV: set AIS_FILE = 'data/ais.csv'
%        CSV columns: Lat, Lon, MMSI (with header)
%     d. Run: >> MAIN_PGA_CFAR
%
%  OUTPUT
%  ------
%  All results are saved to ./results/ and then zipped to:
%    PGA_CFAR_Results.zip   (in the project root)
%
%  Results include:
%    Fig01_Pauli_RGB.png             — Pauli false-color composite
%    Fig02_Polarimetric_Features.png — H, A, alpha, Span, CPD maps
%    Fig03_PSSI_Maps.png             — PSSI components + smoothed PSSI
%    Fig04_PSSI_Gradient_Direction   — Wave propagation direction map
%    Fig05_Adaptive_Parameters.png   — Adaptive k and rho maps
%    Fig06_Detection_Comparison.png  — All 5 detectors side-by-side
%    Fig07_ROC_Curves.png            — PGA-CFAR vs CA-CFAR ROC
%    Fig08_Final_Detection_Map.png   — Final annotated ship positions
%    Fig09_HAlpha_Plane.png          — H/alpha plane scatter
%    Fig10_PSSI_Histogram.png        — PSSI distribution by class
%    Fig11_Ellipse_Kernels.png       — Kernel geometry illustration
%    polarimetric_features.mat       — All feature arrays
%    detection_results.mat           — Detection maps and ship list
%    ablation_results.mat            — Baseline detector maps
%    roc_data.mat                    — ROC sweep data
%    detection_summary.csv           — Ship centroids and areas
%    pssi_statistics.csv             — PSSI statistics by region
%
%  PROJECT STRUCTURE
%  -----------------
%  PGA_CFAR_Project/
%  ├── MAIN_PGA_CFAR.m               Main script (run this)
%  ├── README.m                      This file
%  ├── data/                         Put SAR data here
%  ├── results/                      Auto-created output folder
%  └── src/
%      ├── module1_polsar_reader.m   Data ingestion (real + synthetic)
%      ├── module2_build_coherency.m Coherency matrix T construction
%      ├── module2_lee_filter_T.m    Polarimetric Lee speckle filter
%      ├── module3_feature_engine.m  H, A, alpha, CPD, Bragg ratio
%      ├── module4_pssi_mapper.m     PSSI + gradient (Novel)
%      ├── module5_pga_cfar.m        PGA-CFAR detector (Novel)
%      ├── module6_post_processor.m  Morphological post-processing
%      ├── run_ca_cfar.m             CA-CFAR baseline
%      ├── run_os_cfar.m             OS-CFAR baseline
%      ├── run_h_cfar.m              H-CFAR baseline
%      ├── run_pga_cfar_threshold_only.m  PGA-CFAR thresh-only baseline
%      ├── evaluate_all.m            Pd/Pfa/FoM metrics
%      ├── compute_roc_curves.m      ROC sweep
%      ├── load_ais_csv.m            AIS ground truth loader
%      └── save_all_results.m        Figure generation + file saving
%
%  ALGORITHM PARAMETERS (set in MAIN_PGA_CFAR.m)
%  -----------------------------------------------
%  params.kmin        = 2.5   Threshold scale for calm sea (PSSI=0)
%  params.kmax        = 6.0   Threshold scale for rough sea (PSSI=1)
%  params.r0          = 15    Base ellipse semi-axis (pixels)
%  params.beta        = 3.0   Anisotropy amplification
%  params.alpha_thresh= 40    Pre-screen: min alpha (deg)
%  params.H_thresh    = 0.6   Pre-screen: max entropy
%  params.w1          = 0.4   PSSI weight for entropy H
%  params.w2          = 0.3   PSSI weight for CPD deviation
%  params.w3          = 0.3   PSSI weight for Bragg ratio deviation
%
%  CITATIONS
%  ---------
%  If you use or extend this code, please reference the Phase 1 report:
%  "Sea-State Aware Ship Detection in Ocean Regions Using Polarimetric
%   SAR with a Novel Polarimetric Sea-State Index for Adaptive CFAR
%   Kernel Geometry", Group 11, VIT Chennai, 2024.
% ==========================================================================
