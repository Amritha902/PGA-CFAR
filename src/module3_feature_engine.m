function [H_map, A_map, alpha_map, Span_map, pauli_b2, ...
          CPD_map, Bmeas_map, Bpred_map] = ...
    module3_feature_engine(T, SHH, SHV, SVV, inc_angle_map)
% MODULE3_FEATURE_ENGINE
%   Computes all polarimetric features needed for PSSI and PGA-CFAR:
%
%   FROM EIGENDECOMPOSITION of T:
%     H_map     - Cloude-Pottier entropy H in [0,1]
%     A_map     - Anisotropy A in [0,1]
%     alpha_map - Mean alpha angle in degrees [0,90]
%
%   FROM PAULI DECOMPOSITION:
%     Span_map  - Total power (detection statistic z)
%     pauli_b2  - Double-bounce power |b|^2 (ship indicator)
%
%   FROM CHANNEL DATA:
%     CPD_map   - Co-polarization phase difference angle(T12) [degrees]
%     Bmeas_map - Measured HH/VV backscatter ratio
%     Bpred_map - Bragg-predicted HH/VV ratio from incidence angle

fprintf('  Computing Cloude-Pottier eigendecomposition...\n');
[H_map, A_map, alpha_map] = compute_cloude_pottier(T);

fprintf('  Computing Pauli powers and Span...\n');
Span_map  = T.T11 + T.T22 + T.T33;   % |SHH|^2 + 2|SHV|^2 + |SVV|^2
pauli_b2  = T.T22;                    % |b|^2 = double-bounce power

fprintf('  Computing co-polarization phase difference...\n');
CPD_map = compute_cpd(T);

fprintf('  Computing Bragg ratio...\n');
[Bmeas_map, Bpred_map] = compute_bragg_ratio(T, inc_angle_map);

end


% =========================================================================
function [H, A, alpha] = compute_cloude_pottier(T)
% Pixel-wise eigendecomposition of the 3x3 coherency matrix T.
% Vectorised over all pixels using the analytic 3x3 eigenvalue solver.

[nRows, nCols] = size(T.T11);
H     = zeros(nRows, nCols, 'single');
A     = zeros(nRows, nCols, 'single');
alpha = zeros(nRows, nCols, 'single');

% Build full 3x3 T matrix per pixel and eigen-decompose
% For efficiency, process in blocks
block_size = 128;
for r = 1:block_size:nRows
    r_end = min(r + block_size - 1, nRows);
    for c = 1:block_size:nCols
        c_end = min(c + block_size - 1, nCols);

        for ri = r:r_end
            for ci = c:c_end
                % Assemble T_mat
                T_mat = [T.T11(ri,ci),         T.T12(ri,ci),         T.T13(ri,ci); ...
                         conj(T.T12(ri,ci)),    T.T22(ri,ci),         T.T23(ri,ci); ...
                         conj(T.T13(ri,ci)),    conj(T.T23(ri,ci)),   T.T33(ri,ci)];

                T_mat = (T_mat + T_mat') / 2;  % enforce Hermitian
                T_mat = T_mat + 1e-10*eye(3);  % regularise

                [U, D] = eig(T_mat);
                lambdas = real(diag(D));
                lambdas = max(lambdas, 0);
                [lambdas, idx] = sort(lambdas, 'descend');
                U = U(:, idx);

                sum_L = sum(lambdas) + eps;
                p = lambdas / sum_L;

                % Entropy
                H(ri,ci) = -sum(p .* log(p + eps) / log(3));

                % Anisotropy
                A(ri,ci) = (p(2) - p(3)) / (p(2) + p(3) + eps);

                % Mean alpha angle
                alpha_i = real(acos(min(1, abs(U(1,:)))));  % angle of first component
                alpha(ri,ci) = sum(p .* alpha_i') * (180/pi);
            end
        end
    end
    if mod(r,64) == 1
        fprintf('    Eigendecomp: row %d/%d\n', r_end, nRows);
    end
end
end


% =========================================================================
function CPD = compute_cpd(T)
% Co-polarization phase difference in degrees.
% CPD = angle(<SHH * conj(SVV)>) = angle(T12) in Pauli basis
% In Pauli coherency: T12 = <(SHH+SVV)(SHH-SVV)*> / 2
% The simpler estimate uses the cross-coherence of T elements.
% In the Pauli basis: T12 real/imag encodes the CPD.
% CPD_HH_VV = angle(T11 - T22 + 2*conj(T12)) is exact derivation.
% Practical approach: use angle of <SHH * conj(SVV)> element from C matrix.

% From coherency elements, recover <SHH * SVV*>:
%   C_11 = <|SHH|^2>   C_13 = <SHH * SVV*>   C_33 = <|SVV|^2>
% T12 in Pauli = (<a b*>) where a=(SHH+SVV)/sqrt(2), b=(SHH-SVV)/sqrt(2)
% <SHH * SVV*> = (<aa*> - <bb*>) + j*2*imag(<ab*>)
%              = (T11 - T22)/2 + ... (involves T12)
% Simplest physically correct estimate:
C_SHH_SVV_conj = (T.T11 - T.T22) / 2 - real(T.T12) + 1j*imag(T.T12);
CPD = angle(C_SHH_SVV_conj) * (180/pi);   % degrees
end


% =========================================================================
function [Bmeas, Bpred] = compute_bragg_ratio(T, inc_angle_deg)
% Measured HH/VV ratio and Bragg-predicted ratio.

% Measured ratio from coherency diagonal
sigma_HH = T.T11 + T.T22 + 2*real(T.T12);   % |SHH|^2 from Pauli
sigma_VV = T.T11 + T.T22 - 2*real(T.T12);   % |SVV|^2 from Pauli
Bmeas = (sigma_HH + eps) ./ (sigma_VV + eps);

% Bragg-predicted ratio
theta_r = deg2rad(inc_angle_deg);
eps_r = 70 - 1j*40;   % Seawater at C-band (~5.4 GHz)

cos_t  = cos(theta_r);
sin_t  = sin(theta_r);
eps_cos = sqrt(eps_r - sin_t.^2);

R_HH = (cos_t - eps_cos) ./ (cos_t + eps_cos);
R_VV = (eps_r .* cos_t - eps_cos) ./ (eps_r .* cos_t + eps_cos);

Bpred = abs(R_HH ./ R_VV).^2;
end
