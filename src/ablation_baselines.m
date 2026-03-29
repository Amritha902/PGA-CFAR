% =========================================================================
%  ABLATION BASELINE DETECTORS
%  Implements CA-CFAR, OS-CFAR, H-CFAR, and PGA-CFAR (threshold only)
%  for comparison with full PGA-CFAR (Sec. VIII-B)
% =========================================================================

function [det_map, thr_map] = run_ca_cfar(Span_map, params)
% CA-CFAR: fixed square window, fixed k (kmax as conservative baseline)
k_fixed = (params.kmin + params.kmax) / 2;
r0      = params.r0;
win     = 2*r0 + 1;
guard   = floor(r0/2);

fprintf('    CA-CFAR: k=%.2f, window=%dx%d, guard=%d\n', k_fixed, win, win, guard);

[nR, nC] = size(Span_map);
thr_map  = zeros(nR, nC, 'single');
det_map  = false(nR, nC);

pad = r0 + 1;
Sp  = padarray(Span_map, [pad pad], 'replicate');

for ri = 1:nR
    for ci = 1:nC
        r_p = ri + pad; c_p = ci + pad;
        % Full window minus guard minus CUT
        mask_outer = false(2*r0+1+2, 2*r0+1+2);
        mask_outer(:) = true;
        g = guard;
        mask_outer(pad-g:pad+g, pad-g:pad+g) = false;
        % Extract region
        patch = Sp(r_p-r0-1:r_p+r0+1, c_p-r0-1:c_p+r0+1);
        bg    = patch(mask_outer);
        mu    = mean(bg);
        sig   = std(bg);
        thr   = mu + k_fixed * sig;
        thr_map(ri,ci) = thr;
        det_map(ri,ci) = Span_map(ri,ci) > thr;
    end
end
end


function [det_map, thr_map] = run_os_cfar(Span_map, params)
% OS-CFAR: ordered-statistics, 75th percentile, fixed window
k_fixed  = (params.kmin + params.kmax) / 2;
r0       = params.r0;
pct_rank = 0.75;   % kth order statistic rank

fprintf('    OS-CFAR: k=%.2f, percentile=%.0f%%\n', k_fixed, 100*pct_rank);

[nR, nC] = size(Span_map);
thr_map  = zeros(nR, nC, 'single');
det_map  = false(nR, nC);

pad = r0 + 1;
Sp  = padarray(Span_map, [pad pad], 'replicate');
guard = floor(r0/2);

for ri = 1:nR
    for ci = 1:nC
        r_p = ri + pad; c_p = ci + pad;
        patch = Sp(r_p-r0:r_p+r0, c_p-r0:c_p+r0);
        patch(floor(r0/2):floor(r0/2)+guard+1, ...
              floor(r0/2):floor(r0/2)+guard+1) = NaN;
        bg  = patch(~isnan(patch));
        bg  = sort(bg);
        n   = length(bg);
        kth = bg(min(n, round(pct_rank * n)));
        thr = k_fixed * kth;
        thr_map(ri,ci) = thr;
        det_map(ri,ci) = Span_map(ri,ci) > thr;
    end
end
end


function [det_map, thr_map] = run_h_cfar(Span_map, H_map, params)
% H-CFAR: entropy-adaptive k, fixed square window
% k(m,n) = kmin + (kmax-kmin)*H(m,n)
k_map = params.kmin + (params.kmax - params.kmin) .* H_map;
r0    = params.r0;

fprintf('    H-CFAR: entropy-adaptive k, fixed window\n');

[nR, nC] = size(Span_map);
thr_map  = zeros(nR, nC, 'single');
det_map  = false(nR, nC);
pad      = r0 + 1;
Sp       = padarray(Span_map, [pad pad], 'replicate');
guard    = floor(r0/2);

for ri = 1:nR
    for ci = 1:nC
        r_p = ri + pad; c_p = ci + pad;
        patch = Sp(r_p-r0:r_p+r0, c_p-r0:c_p+r0);
        cx = floor(r0)+1; cy = cx;
        patch(cx-guard:cx+guard, cy-guard:cy+guard) = NaN;
        bg  = patch(~isnan(patch));
        mu  = mean(bg); sig = std(bg);
        thr = mu + k_map(ri,ci) * sig;
        thr_map(ri,ci) = thr;
        det_map(ri,ci) = Span_map(ri,ci) > thr;
    end
end
end


function [det_map, thr_map] = run_pga_cfar_threshold_only(Span_map, PSSI_map, ...
                                    alpha_map, H_map, params)
% PGA-CFAR (threshold only): PSSI-adaptive k + pre-screen, but FIXED window
fprintf('    PGA-CFAR (threshold only): PSSI k, fixed window\n');

k_map   = params.kmin + (params.kmax - params.kmin) .* PSSI_map;
prescreen = (alpha_map > params.alpha_thresh) & (H_map < params.H_thresh);

[nR, nC] = size(Span_map);
thr_map  = zeros(nR, nC, 'single');
det_map  = false(nR, nC);
r0       = params.r0;
pad      = r0 + 1;
Sp       = padarray(Span_map, [pad pad], 'replicate');
guard    = floor(r0/2);

for ri = 1:nR
    for ci = 1:nC
        if ~prescreen(ri,ci), continue; end
        r_p = ri + pad; c_p = ci + pad;
        patch = Sp(r_p-r0:r_p+r0, c_p-r0:c_p+r0);
        cx = floor(r0)+1; cy = cx;
        patch(cx-guard:cx+guard, cy-guard:cy+guard) = NaN;
        bg  = patch(~isnan(patch));
        mu  = mean(bg); sig = std(bg);
        thr = mu + k_map(ri,ci) * sig;
        thr_map(ri,ci) = thr;
        det_map(ri,ci) = Span_map(ri,ci) > thr;
    end
end
end
