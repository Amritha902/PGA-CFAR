function [PSSI, grad_mag, grad_dir] = ...
    module4_pssi_mapper(H_map, CPD_map, Bmeas_map, Bpred_map, params)
% MODULE4_PSSI_MAPPER
%   Computes the Polarimetric Sea-State Index (PSSI) per Eq.(14)-(19):
%
%   PSSI = w1*H + w2*Phi + w3*DeltaB
%
%   where:
%     Phi     = |CPD_HH-VV| / pi  (Eq.15)
%     DeltaB  = min(1, |Bmeas - Bpred| / Bpred)  (Eq.16)
%
%   Then applies Gaussian smoothing and Sobel gradient.
%
%   Outputs:
%     PSSI      - Per-pixel sea-state index in [0,1]
%     grad_mag  - ||nabla PSSI|| (Eq.19)
%     grad_dir  - theta_nabla   (Eq.18) in radians

w1 = params.w1;
w2 = params.w2;
w3 = params.w3;

%% Component 1: Entropy H (already in [0,1])
C1 = double(H_map);
C1 = max(0, min(1, C1));

%% Component 2: Normalized CPD deviation (Eq.15)
CPD_rad = CPD_map * (pi/180);   % to radians
Phi = abs(CPD_rad) / pi;         % [0,1]; Phi->0 Bragg, Phi->1 ship/rough
Phi = max(0, min(1, Phi));

%% Component 3: Bragg ratio deviation (Eq.16)
DeltaB = abs(Bmeas_map - Bpred_map) ./ (Bpred_map + eps);
DeltaB = min(1, DeltaB);
DeltaB = max(0, DeltaB);

%% Weighted combination (Eq.14)
PSSI_raw = w1 * C1 + w2 * Phi + w3 * DeltaB;
PSSI_raw = max(0, min(1, PSSI_raw));

%% Gaussian smoothing to reduce speckle in PSSI map (Sec. VI-D)
sigma_s  = params.pssi_smooth_sigma;
win_s    = params.pssi_smooth_win;
g_kern   = fspecial('gaussian', [win_s win_s], sigma_s);
PSSI     = imfilter(PSSI_raw, g_kern, 'replicate', 'same');

%% Sobel gradient (Eq.17-19)
Kx = [-1 0 1; -2 0 2; -1 0 1];   % Sobel horizontal (range direction)
Ky = [-1 -2 -1; 0 0 0; 1 2 1];   % Sobel vertical (azimuth direction)

dPSSI_dc = imfilter(PSSI, Kx, 'replicate', 'same');  % d/dc (column = range)
dPSSI_dr = imfilter(PSSI, Ky, 'replicate', 'same');  % d/dr (row = azimuth)

grad_mag = sqrt(dPSSI_dc.^2 + dPSSI_dr.^2);          % Eq.19
grad_dir = atan2(dPSSI_dr, dPSSI_dc);                 % Eq.18 (radians)

fprintf('  PSSI components:\n');
fprintf('    C1 (H):      mean=%.3f  std=%.3f\n', mean(C1(:)), std(C1(:)));
fprintf('    C2 (Phi):    mean=%.3f  std=%.3f\n', mean(Phi(:)), std(Phi(:)));
fprintf('    C3 (DeltaB): mean=%.3f  std=%.3f\n', mean(DeltaB(:)), std(DeltaB(:)));
fprintf('  PSSI gradient: max=%.4f  mean=%.4f\n', ...
    max(grad_mag(:)), mean(grad_mag(:)));
end
