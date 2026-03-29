function roc_data = compute_roc_curves(Span_map, PSSI_map, grad_mag, ...
    grad_dir, alpha_map, H_map, params, ship_map_ref)
%COMPUTE_ROC_CURVES Sweeps kmin and computes Pd/Pfa for PGA-CFAR vs CA-CFAR.

kmin_vals = linspace(1.0, 8.0, 15);
n_k = length(kmin_vals);

Pd_pga = zeros(n_k,1); Pfa_pga = zeros(n_k,1);
Pd_ca  = zeros(n_k,1); Pfa_ca  = zeros(n_k,1);

ref_mask = imdilate(logical(ship_map_ref), strel('disk',3));

fprintf('  ROC sweep: %d kmin values [%.1f ... %.1f]\n', ...
    n_k, kmin_vals(1), kmin_vals(end));

for ki = 1:n_k
    p = params;
    p.kmin = kmin_vals(ki);
    p.kmax = kmin_vals(ki) + (params.kmax - params.kmin);

    % PGA-CFAR
    [det_p,~,~,~] = module5_pga_cfar(Span_map, PSSI_map, grad_mag, ...
                        grad_dir, alpha_map, H_map, p);
    [det_p2,~,~,~] = module6_post_processor(det_p, 1);
    det_p2 = logical(det_p2);
    Pd_pga(ki)  = sum(det_p2(:) & ref_mask(:)) / (sum(ref_mask(:)) + eps);
    Pfa_pga(ki) = sum(det_p2(:) & ~ref_mask(:)) / (sum(~ref_mask(:)) + eps);

    % CA-CFAR
    p_ca       = params;
    p_ca.kmin  = kmin_vals(ki);
    p_ca.kmax  = kmin_vals(ki);
    [det_ca,~] = run_ca_cfar(Span_map, p_ca);
    det_ca = logical(det_ca);
    Pd_ca(ki)  = sum(det_ca(:) & ref_mask(:)) / (sum(ref_mask(:)) + eps);
    Pfa_ca(ki) = sum(det_ca(:) & ~ref_mask(:)) / (sum(~ref_mask(:)) + eps);

    fprintf('    k=%.2f  |  PGA: Pd=%.3f Pfa=%.5f  |  CA: Pd=%.3f Pfa=%.5f\n', ...
        kmin_vals(ki), Pd_pga(ki), Pfa_pga(ki), Pd_ca(ki), Pfa_ca(ki));
end

roc_data = struct('kmin_vals', kmin_vals, ...
                  'Pd_pga', Pd_pga, 'Pfa_pga', Pfa_pga, ...
                  'Pd_ca',  Pd_ca,  'Pfa_ca',  Pfa_ca);
end
