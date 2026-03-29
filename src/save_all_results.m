function save_all_results(RESULTS_DIR, ...
    SHH, SHV, SVV, inc_angle_map, ...
    H_map, A_map, alpha_map, Span_map, pauli_b2, CPD_map, ...
    PSSI_map, PSSI_grad_mag, PSSI_grad_dir, ...
    detection_map, ship_map, ship_centroids, ship_areas, ...
    k_map, rho_map, ...
    det_ca, det_os, det_h, det_pga_thresh, ...
    metrics, roc_data, params)
% SAVE_ALL_RESULTS
%   Generates and saves all figures and MAT data files to RESULTS_DIR.

fprintf('  Saving MAT data...\n');
save(fullfile(RESULTS_DIR, 'polarimetric_features.mat'), ...
    'H_map','A_map','alpha_map','Span_map','pauli_b2','CPD_map', ...
    'PSSI_map','PSSI_grad_mag','PSSI_grad_dir','inc_angle_map');
save(fullfile(RESULTS_DIR, 'detection_results.mat'), ...
    'detection_map','ship_map','ship_centroids','ship_areas', ...
    'k_map','rho_map');
save(fullfile(RESULTS_DIR, 'ablation_results.mat'), ...
    'det_ca','det_os','det_h','det_pga_thresh');
save(fullfile(RESULTS_DIR, 'roc_data.mat'), 'roc_data');
if ~isempty(metrics)
    save(fullfile(RESULTS_DIR, 'metrics.mat'), 'metrics');
end
save(fullfile(RESULTS_DIR, 'params.mat'), 'params');

%% ---- FIGURE 1: Pauli RGB False-Color Composite ----
fprintf('  Figure 1: Pauli RGB composite...\n');
fig1 = figure('Visible','off','Position',[0 0 800 800]);
% {R,G,B} = {|b|^2, |c|^2, |a|^2} = {double-bounce, vol, surface}
T11 = real(SHH.*conj(SHH));   % |SHH|^2 ≈ surface+double
T22 = pauli_b2;
T33 = real(SVV.*conj(SVV));
pauli_R = db_scale(T22);
pauli_G = db_scale(real(SHV.*conj(SHV)));
pauli_B = db_scale(real((SHH+SVV).*conj(SHH+SVV))/2);
rgb_pauli = cat(3, pauli_R, pauli_G, pauli_B);
imshow(rgb_pauli);
title('Pauli False-Color Composite (R=Double-bounce, G=Volume, B=Surface)', ...
    'FontSize',11,'FontWeight','bold');
saveas(fig1, fullfile(RESULTS_DIR, 'Fig01_Pauli_RGB.png'));
close(fig1);

%% ---- FIGURE 2: Polarimetric Feature Maps ----
fprintf('  Figure 2: Polarimetric feature maps...\n');
fig2 = figure('Visible','off','Position',[0 0 1200 800]);
subplot(2,3,1); imagesc(H_map,[0 1]);      colorbar; title('Entropy H'); colormap(gca,'jet');
subplot(2,3,2); imagesc(A_map,[0 1]);      colorbar; title('Anisotropy A'); colormap(gca,'jet');
subplot(2,3,3); imagesc(alpha_map,[0 90]); colorbar; title('Mean Alpha (deg)'); colormap(gca,'jet');
subplot(2,3,4); imagesc(10*log10(Span_map+eps),[-25 5]); colorbar; title('Span (dB)'); colormap(gca,'gray');
subplot(2,3,5); imagesc(CPD_map,[-180 180]); colorbar; title('CPD \phi_{HH-VV} (deg)'); colormap(gca,'hsv');
subplot(2,3,6); imagesc(inc_angle_map);   colorbar; title('Incidence Angle (deg)'); colormap(gca,'jet');
sgtitle('Polarimetric Feature Maps — Stage 3', 'FontSize',13,'FontWeight','bold');
saveas(fig2, fullfile(RESULTS_DIR, 'Fig02_Polarimetric_Features.png'));
close(fig2);

%% ---- FIGURE 3: PSSI Map and Components ----
fprintf('  Figure 3: PSSI map...\n');
fig3 = figure('Visible','off','Position',[0 0 1200 500]);

% Recompute components for display
C1   = max(0,min(1,double(H_map)));
CPD_rad = CPD_map * (pi/180);
Phi  = max(0,min(1,abs(CPD_rad)/pi));
DeltaB = min(1,max(0, abs(real(ones(size(PSSI_map)))) ));  % approximate

subplot(1,4,1); imagesc(C1,[0 1]);     colorbar; title('C1: Entropy H'); colormap(gca,'hot');
subplot(1,4,2); imagesc(Phi,[0 1]);    colorbar; title('C2: \Phi (CPD dev.)'); colormap(gca,'hot');
subplot(1,4,3); imagesc(PSSI_map,[0 1]); colorbar; title('PSSI (smoothed)'); colormap(gca,'jet');
subplot(1,4,4); imagesc(PSSI_grad_mag); colorbar; title('||\nabla PSSI||'); colormap(gca,'parula');
sgtitle('Polarimetric Sea-State Index (PSSI) — Stage 4', 'FontSize',13,'FontWeight','bold');
saveas(fig3, fullfile(RESULTS_DIR, 'Fig03_PSSI_Maps.png'));
close(fig3);

%% ---- FIGURE 4: PSSI Gradient Direction Map ----
fprintf('  Figure 4: PSSI gradient direction (wave orientation)...\n');
fig4 = figure('Visible','off','Position',[0 0 700 700]);
imagesc(PSSI_grad_dir * 180/pi, [-180 180]);
colorbar; colormap(gca,'hsv');
hold on;
% Overlay gradient arrows (subsampled)
step = 20;
[Y,X] = meshgrid(1:step:size(PSSI_grad_dir,2), 1:step:size(PSSI_grad_dir,1));
U = cos(PSSI_grad_dir(1:step:end, 1:step:end));
V = sin(PSSI_grad_dir(1:step:end, 1:step:end));
quiver(Y, X, U, V, 0.8, 'w', 'LineWidth', 0.5);
title('PSSI Gradient Direction \theta_\nabla (wave propagation)', ...
    'FontSize',12,'FontWeight','bold');
xlabel('Column (Range)'); ylabel('Row (Azimuth)');
saveas(fig4, fullfile(RESULTS_DIR, 'Fig04_PSSI_Gradient_Direction.png'));
close(fig4);

%% ---- FIGURE 5: Adaptive Parameter Maps ----
fprintf('  Figure 5: Adaptive k and rho maps...\n');
fig5 = figure('Visible','off','Position',[0 0 900 400]);
subplot(1,2,1); imagesc(k_map, [params.kmin params.kmax]); colorbar;
title(sprintf('Adaptive k map [%.1f, %.1f]', params.kmin, params.kmax));
colormap(gca,'parula');
subplot(1,2,2); imagesc(rho_map, [1, 1+params.beta]); colorbar;
title(sprintf('Anisotropy \\rho [1, %.1f]', 1+params.beta));
colormap(gca,'hot');
sgtitle('PGA-CFAR Adaptive Parameters — Stage 5', 'FontSize',13,'FontWeight','bold');
saveas(fig5, fullfile(RESULTS_DIR, 'Fig05_Adaptive_Parameters.png'));
close(fig5);

%% ---- FIGURE 6: Detection Results Comparison ----
fprintf('  Figure 6: Detection comparison...\n');
fig6 = figure('Visible','off','Position',[0 0 1400 900]);
Span_dB = 10*log10(Span_map + eps);
clim    = [prctile(Span_dB(:),2), prctile(Span_dB(:),98)];

method_names = {'CA-CFAR', 'OS-CFAR', 'H-CFAR', ...
                'PGA-CFAR\n(Thresh Only)', 'Full PGA-CFAR'};
det_list = {det_ca, det_os, det_h, det_pga_thresh, ship_map};

for i = 1:5
    subplot(2,3,i);
    imagesc(Span_dB, clim); colormap(gca,'gray'); hold on;
    % Mark detections
    [dr, dc] = find(det_list{i});
    if ~isempty(dr)
        plot(dc, dr, 'r.', 'MarkerSize', 3);
    end
    title(method_names{i}, 'FontSize',10,'FontWeight','bold');
    xlabel('Range'); ylabel('Azimuth');
end
subplot(2,3,6);
imagesc(Span_dB, clim); colormap(gca,'gray'); hold on;
if ~isempty(ship_centroids)
    plot(ship_centroids(:,1), ship_centroids(:,2), 'r+', 'MarkerSize',12, 'LineWidth',2);
    plot(ship_centroids(:,1), ship_centroids(:,2), 'yo', 'MarkerSize',15, 'LineWidth',1.5);
end
title('Full PGA-CFAR — Final Ships', 'FontSize',10,'FontWeight','bold');
sgtitle('Ablation Study: Detection Map Comparison', 'FontSize',14,'FontWeight','bold');
saveas(fig6, fullfile(RESULTS_DIR, 'Fig06_Detection_Comparison.png'));
close(fig6);

%% ---- FIGURE 7: ROC Curves ----
fprintf('  Figure 7: ROC curves...\n');
if ~isempty(roc_data)
    fig7 = figure('Visible','off','Position',[0 0 700 600]);
    plot(roc_data.Pfa_ca,  roc_data.Pd_ca,  'b-o', 'LineWidth',2, ...
        'MarkerFaceColor','b', 'DisplayName','CA-CFAR'); hold on;
    plot(roc_data.Pfa_pga, roc_data.Pd_pga, 'r-s', 'LineWidth',2, ...
        'MarkerFaceColor','r', 'DisplayName','Full PGA-CFAR');
    xlabel('False Alarm Rate P_{fa}', 'FontSize',12);
    ylabel('Probability of Detection P_d', 'FontSize',12);
    title('ROC Curve: PGA-CFAR vs CA-CFAR', 'FontSize',13,'FontWeight','bold');
    legend('Location','southeast','FontSize',11);
    grid on; xlim([0 0.05]); ylim([0 1]);
    text(0.03, 0.3, sprintf('PGA-CFAR AUC \\approx %.3f', ...
        trapz(roc_data.Pfa_pga, roc_data.Pd_pga)), 'FontSize',10, 'Color','r');
    saveas(fig7, fullfile(RESULTS_DIR, 'Fig07_ROC_Curves.png'));
    close(fig7);
end

%% ---- FIGURE 8: Final Detection Overlay ----
fprintf('  Figure 8: Final detection overlay...\n');
fig8 = figure('Visible','off','Position',[0 0 900 900]);
imagesc(Span_dB, clim); colormap('gray'); hold on; axis image;
[sr, sc] = find(ship_map);
if ~isempty(sr)
    plot(sc, sr, 'r.', 'MarkerSize', 2);
end
if ~isempty(ship_centroids)
    for i = 1:size(ship_centroids,1)
        plot(ship_centroids(i,1), ship_centroids(i,2), ...
            'r+', 'MarkerSize',14, 'LineWidth',2);
        text(ship_centroids(i,1)+3, ship_centroids(i,2)-3, ...
            sprintf('S%d', i), 'Color','yellow', 'FontSize',9, 'FontWeight','bold');
    end
end
colorbar;
title(sprintf('PGA-CFAR Final Detection Map — %d Ships', size(ship_centroids,1)), ...
    'FontSize',13,'FontWeight','bold');
xlabel('Column (Range)'); ylabel('Row (Azimuth)');
saveas(fig8, fullfile(RESULTS_DIR, 'Fig08_Final_Detection_Map.png'));
close(fig8);

%% ---- FIGURE 9: H/alpha plane ----
fprintf('  Figure 9: H/alpha plane...\n');
fig9 = figure('Visible','off','Position',[0 0 700 600]);
H_flat     = H_map(:);
alpha_flat = alpha_map(:);
ship_flat  = ship_map(:);
scatter(H_flat(~ship_flat), alpha_flat(~ship_flat), 1, [0.5 0.7 1], ...
    'filled', 'DisplayName','Ocean'); hold on;
scatter(H_flat(ship_flat), alpha_flat(ship_flat), 10, 'r', ...
    'filled', 'DisplayName','Detected Ships');
xlabel('Entropy H', 'FontSize',12);
ylabel('Mean Alpha \alpha (deg)', 'FontSize',12);
title('H/\alpha Plane — Ocean vs Detected Ships', 'FontSize',13,'FontWeight','bold');
legend('Location','northwest','FontSize',11);
xlim([0 1]); ylim([0 90]);
grid on;
% Draw Cloude-Pottier zone boundaries (approximate)
xline(0.5,'--','Color',[0.4 0.4 0.4],'LineWidth',1.5);
xline(0.9,'--','Color',[0.4 0.4 0.4],'LineWidth',1.5);
yline(40,':','Color',[0.4 0.4 0.4],'LineWidth',1.5);
saveas(fig9, fullfile(RESULTS_DIR, 'Fig09_HAlpha_Plane.png'));
close(fig9);

%% ---- FIGURE 10: PSSI histogram by region ----
fprintf('  Figure 10: PSSI histogram...\n');
fig10 = figure('Visible','off','Position',[0 0 700 500]);
pssi_ship = PSSI_map(ship_map);
pssi_ocean = PSSI_map(~ship_map);
histogram(pssi_ocean, 50, 'Normalization','probability', ...
    'FaceColor',[0.3 0.6 1], 'EdgeColor','none', 'DisplayName','Ocean'); hold on;
histogram(pssi_ship,  50, 'Normalization','probability', ...
    'FaceColor',[1 0.3 0.3], 'EdgeColor','none', 'DisplayName','Detected Ships');
xlabel('PSSI Value', 'FontSize',12);
ylabel('Normalized Frequency', 'FontSize',12);
title('PSSI Distribution: Ocean vs Detected Ships', 'FontSize',13,'FontWeight','bold');
legend('FontSize',11);
grid on;
saveas(fig10, fullfile(RESULTS_DIR, 'Fig10_PSSI_Histogram.png'));
close(fig10);

%% ---- FIGURE 11: Ellipse kernel illustration ----
fprintf('  Figure 11: Ellipse kernel illustration...\n');
fig11 = figure('Visible','off','Position',[0 0 700 700]);
theta_examples = [0, pi/4, pi/2, 3*pi/4];
rho_examples   = [1.0, 2.0, 3.0, 4.0];
r0_illus = 15;

for i = 1:4
    subplot(2,2,i);
    rho = rho_examples(i);
    th  = theta_examples(i);
    a   = r0_illus * rho; b = r0_illus / rho;
    ag  = a/2; bg = b/2;
    t   = linspace(0, 2*pi, 200);
    % Outer ellipse (rotated)
    xo = a*cos(t); yo = b*sin(t);
    xo_r = xo*cos(th) - yo*sin(th);
    yo_r = xo*sin(th) + yo*cos(th);
    % Inner guard
    xi = ag*cos(t); yi = bg*sin(t);
    xi_r = xi*cos(th) - yi*sin(th);
    yi_r = xi*sin(th) + yi*cos(th);
    plot(xo_r, yo_r, 'b-', 'LineWidth',2); hold on;
    plot(xi_r, yi_r, 'r-', 'LineWidth',2);
    plot(0, 0, 'k+', 'MarkerSize',12, 'LineWidth',2);
    % Wave direction arrow
    arrow_len = r0_illus * 1.5;
    quiver(0,0,arrow_len*cos(th),arrow_len*sin(th),0,'g','LineWidth',2);
    axis equal; grid on;
    title(sprintf('\\rho=%.1f, \\theta=%.0f°', rho, th*180/pi));
    xlabel('Range'); ylabel('Azimuth');
    xlim([-60 60]); ylim([-60 60]);
end
sgtitle('Adaptive Elliptical Kernel Geometry (Blue=BG, Red=Guard)', ...
    'FontSize',12,'FontWeight','bold');
saveas(fig11, fullfile(RESULTS_DIR, 'Fig11_Ellipse_Kernels.png'));
close(fig11);

%% ---- Save numerical summary CSV ----
fprintf('  Saving numerical summary...\n');
fid = fopen(fullfile(RESULTS_DIR, 'detection_summary.csv'), 'w');
fprintf(fid, 'Ship_ID,Centroid_Col_px,Centroid_Row_px,Area_px\n');
for i = 1:size(ship_centroids,1)
    fprintf(fid, '%d,%.1f,%.1f,%d\n', i, ...
        ship_centroids(i,1), ship_centroids(i,2), ship_areas(i));
end
fclose(fid);

%% ---- Save PSSI statistics CSV ----
fid2 = fopen(fullfile(RESULTS_DIR, 'pssi_statistics.csv'), 'w');
fprintf(fid2, 'Region,PSSI_mean,PSSI_std,PSSI_min,PSSI_max\n');
fprintf(fid2, 'Full Scene,%.4f,%.4f,%.4f,%.4f\n', ...
    mean(PSSI_map(:)), std(PSSI_map(:)), min(PSSI_map(:)), max(PSSI_map(:)));
if any(ship_map(:))
    fprintf(fid2, 'Ship Pixels,%.4f,%.4f,%.4f,%.4f\n', ...
        mean(PSSI_map(ship_map)), std(PSSI_map(ship_map)), ...
        min(PSSI_map(ship_map)), max(PSSI_map(ship_map)));
    fprintf(fid2, 'Ocean Pixels,%.4f,%.4f,%.4f,%.4f\n', ...
        mean(PSSI_map(~ship_map)), std(PSSI_map(~ship_map)), ...
        min(PSSI_map(~ship_map)), max(PSSI_map(~ship_map)));
end
fclose(fid2);

fprintf('  All results saved to: %s\n', RESULTS_DIR);
fprintf('  Files generated:\n');
d = dir(fullfile(RESULTS_DIR,'*'));
for i = 1:length(d)
    if ~d(i).isdir
        fprintf('    %s (%.1f KB)\n', d(i).name, d(i).bytes/1024);
    end
end
end


% =========================================================================
%  Helper: log-scale stretch for display
% =========================================================================
function out = db_scale(x)
    x_db = 10*log10(x + eps);
    lo = prctile(x_db(:), 2);
    hi = prctile(x_db(:), 98);
    out = (x_db - lo) / (hi - lo + eps);
    out = max(0, min(1, out));
end
