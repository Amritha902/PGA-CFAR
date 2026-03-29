function metrics = evaluate_all(ship_map_pga, det_ca, det_os, det_h, ...
                                det_pga_thresh, ais_data)
%EVALUATE_ALL Computes Pd, Pfa, FoM for each detector vs AIS ground truth.

methods   = {'CA_CFAR','OS_CFAR','H_CFAR','PGA_thresh_only','Full_PGA_CFAR'};
det_maps  = {det_ca, det_os, det_h, det_pga_thresh, ship_map_pga};

fprintf('\n  %-26s %8s %10s %10s\n','Method','Pd','Pfa','FoM');
fprintf('  %s\n', repmat('-',1,58));

metrics = struct();
ship_radius = 3;
n_ships = length(ais_data.row_px);
[nR,nC] = size(ship_map_pga);

for mi = 1:length(methods)
    det = logical(det_maps{mi});
    n_det = 0;
    for s = 1:n_ships
        r0 = ais_data.row_px(s); c0 = ais_data.col_px(s);
        rr = max(1,r0-ship_radius):min(nR,r0+ship_radius);
        cc = max(1,c0-ship_radius):min(nC,c0+ship_radius);
        if any(any(det(rr,cc))), n_det = n_det+1; end
    end
    Pd  = n_det / (n_ships+eps);
    n_fp = max(0, sum(det(:)) - n_det*(2*ship_radius+1)^2);
    Pfa = n_fp / (nR*nC);
    FoM = Pd - 10*Pfa;
    metrics.(methods{mi}).Pd  = Pd;
    metrics.(methods{mi}).Pfa = Pfa;
    metrics.(methods{mi}).FoM = FoM;
    fprintf('  %-26s %8.3f %10.5f %10.3f\n', strrep(methods{mi},'_',' '), Pd, Pfa, FoM);
end
end
