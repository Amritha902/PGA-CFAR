% =========================================================================
%  MAIN_PGA_CFAR.m
%  Sea-State Aware Ship Detection using Polarimetric SAR
%  PGA-CFAR Framework — Group 11, VIT Chennai
%
%  USAGE:
%    1. Set DATA_PATH to your quad-pol SAR data folder
%    2. Set DATA_FORMAT to your sensor type
%    3. (Optional) Set AIS_FILE for ground truth evaluation
%    4. Run this script: >> MAIN_PGA_CFAR
%
%  OUTPUT:
%    All figures and result files are saved to ./results/
%    A ZIP archive (PGA_CFAR_Results.zip) is created automatically.
% =========================================================================

clc; clear; close all;
fprintf('=================================================================\n');
fprintf('  PGA-CFAR Ship Detection Framework\n');
fprintf('  Group 11 — VIT Chennai\n');
fprintf('=================================================================\n\n');

%% -------------------------------------------------------------------------
%  USER CONFIGURATION — EDIT THESE PATHS BEFORE RUNNING
% -------------------------------------------------------------------------
% Path to the folder containing your SAR data files
DATA_PATH   = './data';

% Sensor / format selector:
%   'RADARSAT2'  -> reads SHH.bin, SHV.bin, SVV.bin + lutSigma.xml
%   'ALOS2'      -> reads CEOS format (IMG-HH-*, IMG-HV-*, IMG-VV-*)
%   'GAOFEN3'    -> reads GF3 GeoTIFF quad-pol
%   'SYNTHETIC'  -> generates synthetic scene (for code testing only)
DATA_FORMAT = 'SYNTHETIC';

% (Optional) Path to AIS ground-truth CSV: columns [Lat, Lon, MMSI]
% Set to '' to skip evaluation
AIS_FILE    = '';

% Results output folder
RESULTS_DIR = './results';

%% -------------------------------------------------------------------------
%  PGA-CFAR ALGORITHM PARAMETERS (from paper)
% -------------------------------------------------------------------------
params.kmin        = 2.5;    % Minimum CFAR threshold scale (calm sea)
params.kmax        = 6.0;    % Maximum CFAR threshold scale (rough sea)
params.r0          = 15;     % Base ellipse semi-axis radius (pixels)
params.beta        = 3.0;    % Anisotropy amplification factor
params.alpha_thresh= 40;     % Pre-screen: min alpha angle (degrees)
params.H_thresh    = 0.6;    % Pre-screen: max entropy for ship candidates
params.w1          = 0.4;    % PSSI weight: Entropy H
params.w2          = 0.3;    % PSSI weight: CPD deviation
params.w3          = 0.3;    % PSSI weight: Bragg ratio deviation
params.pssi_smooth_sigma = 2;% Gaussian smoothing sigma for PSSI (pixels)
params.pssi_smooth_win   = 11;% Gaussian smoothing window size
params.T_avg_win   = 7;      % Coherency matrix spatial averaging window
params.lee_win     = 5;      % Refined Lee filter window size
params.min_ship_area = 4;    % Minimum ship component area (pixels^2)

%% -------------------------------------------------------------------------
%  SETUP
% -------------------------------------------------------------------------
if ~exist(RESULTS_DIR, 'dir'), mkdir(RESULTS_DIR); end
addpath(genpath('./src'));
rng(42);   % reproducibility for synthetic scenes

fprintf('[SETUP] Results will be saved to: %s\n\n', RESULTS_DIR);

%% =========================================================================
%  MODULE 1: DATA INGESTION
% =========================================================================
fprintf('--- MODULE 1: Loading PolSAR Data ---\n');
[SHH, SHV, SVV, inc_angle_map, geo_info] = ...
    module1_polsar_reader(DATA_PATH, DATA_FORMAT);

[Nrows, Ncols] = size(SHH);
fprintf('  Scene size: %d x %d pixels\n', Nrows, Ncols);
fprintf('  Incidence angle range: %.1f - %.1f deg\n', ...
    min(inc_angle_map(:)), max(inc_angle_map(:)));

%% =========================================================================
%  MODULE 2: PREPROCESSING
% =========================================================================
fprintf('\n--- MODULE 2: Preprocessing ---\n');

% 2a. Radiometric calibration (already done in reader for real data;
%     for synthetic data, SHH/SHV/SVV are already calibrated)
fprintf('  Radiometric calibration: done\n');

% 2b. Build coherency matrix T with 7x7 spatial averaging
fprintf('  Building coherency matrix T (7x7 averaging)...\n');
T = module2_build_coherency(SHH, SHV, SVV, params.T_avg_win);

% 2c. Polarimetric-preserving Lee filter on T
fprintf('  Applying refined Lee speckle filter (%dx%d)...\n', ...
    params.lee_win, params.lee_win);
T = module2_lee_filter_T(T, params.lee_win);

fprintf('  Preprocessing complete.\n');

%% =========================================================================
%  MODULE 3: POLARIMETRIC FEATURE EXTRACTION
% =========================================================================
fprintf('\n--- MODULE 3: Polarimetric Feature Extraction ---\n');

[H_map, A_map, alpha_map, Span_map, pauli_b2, ...
    CPD_map, Bmeas_map, Bpred_map] = ...
    module3_feature_engine(T, SHH, SHV, SVV, inc_angle_map);

fprintf('  H:     mean=%.3f, std=%.3f\n', mean(H_map(:)), std(H_map(:)));
fprintf('  alpha: mean=%.1f deg\n', mean(alpha_map(:)));
fprintf('  Span:  mean=%.2f dB\n', 10*log10(mean(Span_map(:))+eps));

%% =========================================================================
%  MODULE 4: PSSI COMPUTATION
% =========================================================================
fprintf('\n--- MODULE 4: PSSI Computation ---\n');

[PSSI_map, PSSI_grad_mag, PSSI_grad_dir] = ...
    module4_pssi_mapper(H_map, CPD_map, Bmeas_map, Bpred_map, params);

fprintf('  PSSI range: [%.3f, %.3f]\n', min(PSSI_map(:)), max(PSSI_map(:)));
fprintf('  PSSI mean:  %.3f\n', mean(PSSI_map(:)));

%% =========================================================================
%  MODULE 5: PGA-CFAR DETECTION
% =========================================================================
fprintf('\n--- MODULE 5: PGA-CFAR Detection ---\n');

[detection_map, k_map, rho_map, ellipse_params] = ...
    module5_pga_cfar(Span_map, PSSI_map, PSSI_grad_mag, PSSI_grad_dir, ...
                     alpha_map, H_map, params);

n_candidates = sum(detection_map(:));
fprintf('  Raw detections (before post-processing): %d pixels\n', n_candidates);

%% =========================================================================
%  MODULE 6: POST-PROCESSING
% =========================================================================
fprintf('\n--- MODULE 6: Post-Processing ---\n');

[ship_map, ship_centroids, ship_areas, ship_labels] = ...
    module6_post_processor(detection_map, params.min_ship_area);

n_ships = size(ship_centroids, 1);
fprintf('  Ships detected after post-processing: %d\n', n_ships);
if n_ships > 0
    fprintf('  Ship areas (px): min=%d, max=%d, mean=%.1f\n', ...
        min(ship_areas), max(ship_areas), mean(ship_areas));
end

%% =========================================================================
%  MODULE 7: EVALUATION & ABLATION
% =========================================================================
fprintf('\n--- MODULE 7: Evaluation & Ablation Study ---\n');

% Run ablation baselines for comparison
fprintf('  Running CA-CFAR baseline...\n');
[det_ca, ~] = run_ca_cfar(Span_map, params);

fprintf('  Running OS-CFAR baseline...\n');
[det_os, ~] = run_os_cfar(Span_map, params);

fprintf('  Running H-CFAR (entropy-only) baseline...\n');
[det_h, ~] = run_h_cfar(Span_map, H_map, params);

fprintf('  Running PGA-CFAR (threshold only, fixed window)...\n');
[det_pga_thresh, ~] = run_pga_cfar_threshold_only(Span_map, PSSI_map, ...
                        alpha_map, H_map, params);

% Evaluate if AIS available
if ~isempty(AIS_FILE) && exist(AIS_FILE, 'file')
    fprintf('  Loading AIS ground truth: %s\n', AIS_FILE);
    ais_data = load_ais_csv(AIS_FILE, geo_info, Nrows, Ncols);
    metrics = evaluate_all(ship_map, det_ca, det_os, det_h, ...
                           det_pga_thresh, ais_data);
else
    fprintf('  No AIS file provided — skipping quantitative evaluation.\n');
    fprintf('  Set AIS_FILE in the configuration block to enable metrics.\n');
    metrics = [];
end

% Compute ROC curve (over kmin sweep)
fprintf('  Computing ROC curves...\n');
[roc_data] = compute_roc_curves(Span_map, PSSI_map, PSSI_grad_mag, ...
    PSSI_grad_dir, alpha_map, H_map, params, ship_map);

%% =========================================================================
%  SAVE ALL RESULTS
% =========================================================================
fprintf('\n--- Saving Results ---\n');
save_all_results(RESULTS_DIR, ...
    SHH, SHV, SVV, inc_angle_map, ...
    H_map, A_map, alpha_map, Span_map, pauli_b2, CPD_map, ...
    PSSI_map, PSSI_grad_mag, PSSI_grad_dir, ...
    detection_map, ship_map, ship_centroids, ship_areas, ...
    k_map, rho_map, ...
    det_ca, det_os, det_h, det_pga_thresh, ...
    metrics, roc_data, params);

%% =========================================================================
%  CREATE ZIP ARCHIVE
% =========================================================================
fprintf('\n--- Creating ZIP Archive ---\n');
zip_path = fullfile(pwd, 'PGA_CFAR_Results.zip');
zip(zip_path, RESULTS_DIR);
fprintf('  ZIP saved to: %s\n', zip_path);

fprintf('\n=================================================================\n');
fprintf('  ALL DONE. Download PGA_CFAR_Results.zip\n');
fprintf('=================================================================\n');
