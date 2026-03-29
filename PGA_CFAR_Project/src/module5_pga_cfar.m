function [detection_map, k_map, rho_map, ellipse_params] = ...
    module5_pga_cfar(Span_map, PSSI_map, grad_mag, grad_dir, ...
                     alpha_map, H_map, params)
% MODULE5_PGA_CFAR
%   Implements the PGA-CFAR detection algorithm (Algorithm 1).
%   Both the threshold scaling factor k and the CFAR window geometry
%   are adapted pixel-by-pixel from the PSSI map.
%
%   Decision rule (Eq.26):
%     y(m,n) = 1  if z(m,n) > mu_R(m,n) + k(m,n)*sigma_R(m,n)
%
%   Novel contributions:
%     1. Adaptive k from PSSI (Eq.21)
%     2. Adaptive elliptical kernel oriented along wave crests (Eq.23-25)

kmin   = params.kmin;
kmax   = params.kmax;
r0     = params.r0;
beta   = params.beta;
a_thr  = params.alpha_thresh;
H_thr  = params.H_thresh;

[nRows, nCols] = size(Span_map);
detection_map  = false(nRows, nCols);

%% --- Step 1: Adaptive k map (Eq.21) ---
k_map = kmin + (kmax - kmin) .* PSSI_map;

%% --- Step 2: Anisotropy factor rho (Eq.23) ---
% Normalise gradient magnitude to [0,1] for rho computation
grad_norm = grad_mag / (max(grad_mag(:)) + eps);
rho_map   = 1 + beta .* grad_norm;

%% --- Step 3: Ellipse axes ---
a_map = r0 .* rho_map;          % semi-major (along wave crest)
b_map = r0 ./ rho_map;          % semi-minor (wave propagation dir)
% a*b = r0^2 = constant (total background cells approximately fixed)

% Guard ellipse (half size)
ag_map = a_map / 2;
bg_map = b_map / 2;

ellipse_params = struct('a', a_map, 'b', b_map, 'theta', grad_dir, ...
                        'rho', rho_map);

%% --- Step 4: Polarimetric pre-screen (Sec. V-D) ---
% Candidate pixels must satisfy: alpha > alpha_thresh AND H < H_thresh
prescreen_mask = (alpha_map > a_thr) & (H_map < H_thr);
fprintf('  Pre-screen: %.2f%% pixels are candidates\n', ...
    100*mean(prescreen_mask(:)));

%% --- Step 5: Per-pixel CFAR ---
% Precompute a lookup table of ellipse membership for discrete (rho, theta) bins
% to avoid per-pixel mask recomputation (as described in Sec. IX)
fprintf('  Building ellipse LUT...\n');
[lut_rho_vals, lut_theta_vals, lut_cells] = build_ellipse_lut(r0, beta);
fprintf('  Running PGA-CFAR on %d candidate pixels...\n', sum(prescreen_mask(:)));

% Pad Span image to avoid border issues
pad_size = ceil(r0 * (1 + beta)) + 2;
Span_pad = padarray(Span_map, [pad_size pad_size], 'replicate', 'both');

n_detected = 0;
candidates = find(prescreen_mask);
total_cand = length(candidates);

report_interval = max(1, round(total_cand / 10));

for idx = 1:total_cand
    if mod(idx, report_interval) == 0
        fprintf('    Progress: %d/%d (%.0f%%)\n', idx, total_cand, ...
            100*idx/total_cand);
    end

    [ri, ci] = ind2sub([nRows nCols], candidates(idx));

    % Lookup LUT bin for this pixel's rho and theta
    rho_val   = rho_map(ri, ci);
    theta_val = grad_dir(ri, ci);
    bg_cells  = lut_lookup(lut_rho_vals, lut_theta_vals, lut_cells, ...
                           rho_val, theta_val, r0);

    % Offset background cells into padded image coordinates
    bg_r = ri + pad_size + bg_cells(:,1);
    bg_c = ci + pad_size + bg_cells(:,2);

    % Clip to padded image bounds
    valid = bg_r >= 1 & bg_r <= size(Span_pad,1) & ...
            bg_c >= 1 & bg_c <= size(Span_pad,2);
    bg_r = bg_r(valid); bg_c = bg_c(valid);

    if length(bg_r) < 4
        continue;  % Not enough background cells
    end

    bg_vals = Span_pad(sub2ind(size(Span_pad), bg_r, bg_c));
    mu_R    = mean(bg_vals);
    sig_R   = std(bg_vals);

    z   = Span_map(ri, ci);
    k   = k_map(ri, ci);
    thr = mu_R + k * sig_R;

    if z > thr
        detection_map(ri, ci) = true;
        n_detected = n_detected + 1;
    end
end

fprintf('  PGA-CFAR: %d raw detection pixels\n', n_detected);
end


% =========================================================================
%  ELLIPSE LOOKUP TABLE BUILDER
% =========================================================================
function [rho_vals, theta_vals, cell_list] = build_ellipse_lut(r0, beta)
% Discretise (rho, theta) space and precompute background cell offsets.

rho_vals   = linspace(1.0, 1 + beta, 12);    % 12 rho bins
theta_vals = linspace(-pi, pi, 18);           % 18 angle bins (20-degree steps)

cell_list = cell(length(rho_vals), length(theta_vals));

max_r = ceil(r0 * (1 + beta)) + 2;
[dx, dy] = meshgrid(-max_r:max_r, -max_r:max_r);
dx = dx(:); dy = dy(:);

for ri = 1:length(rho_vals)
    rho = rho_vals(ri);
    a   = r0 * rho;
    b   = r0 / rho;
    ag  = a / 2;
    bg  = b / 2;

    for ti = 1:length(theta_vals)
        theta = theta_vals(ti);
        ct = cos(theta); st = sin(theta);

        % Rotated ellipse membership (Eq.24)
        x_rot =  dx*ct + dy*st;
        y_rot = -dx*st + dy*ct;

        outer = (x_rot/a).^2 + (y_rot/b).^2 <= 1;
        inner = (x_rot/ag).^2 + (y_rot/bg).^2 <= 1;

        bg_mask = outer & ~inner;
        cell_list{ri, ti} = [dy(bg_mask), dx(bg_mask)];  % [row_offset, col_offset]
    end
end
end


function cells = lut_lookup(rho_vals, theta_vals, lut, rho, theta, ~)
% Find nearest LUT bin and return precomputed cell offsets.
    [~, ri] = min(abs(rho_vals - rho));
    [~, ti] = min(abs(theta_vals - theta));
    cells   = lut{ri, ti};
end
