function T = module2_build_coherency(SHH, SHV, SVV, win_size)
% MODULE2_BUILD_COHERENCY
%   Constructs the 3x3 Hermitian coherency matrix T(m,n) from complex
%   scattering matrix elements using Pauli decomposition and spatial
%   averaging over a win_size x win_size boxcar window.
%
%   T = < kP * kP^H > where kP = 1/sqrt(2) * [SHH+SVV; SHH-SVV; 2*SHV]
%
%   Input:
%     SHH, SHV, SVV  : Complex scattering matrix channels (NxM)
%     win_size        : Spatial averaging window size (pixels)
%   Output:
%     T               : Coherency matrix struct with fields T11..T33
%                       (each NxM real or complex matrix)

% Pauli basis components
a = (SHH + SVV) / sqrt(2);   % Surface scattering
b = (SHH - SVV) / sqrt(2);   % Double-bounce
c = sqrt(2) * SHV;            % Volume / cross

% All 9 outer products of kP = [a; b; c]
% kP * kP^H elements (use spatial average as ensemble average)
T.T11 = boxcar_avg(a .* conj(a), win_size);  % Real
T.T22 = boxcar_avg(b .* conj(b), win_size);  % Real
T.T33 = boxcar_avg(c .* conj(c), win_size);  % Real
T.T12 = boxcar_avg(a .* conj(b), win_size);  % Complex
T.T13 = boxcar_avg(a .* conj(c), win_size);  % Complex
T.T23 = boxcar_avg(b .* conj(c), win_size);  % Complex
% Lower triangle by Hermitian symmetry
% T21 = conj(T12), T31 = conj(T13), T32 = conj(T23)

T.T11 = real(T.T11);
T.T22 = real(T.T22);
T.T33 = real(T.T33);

end


% =========================================================================
function T = module2_lee_filter_T(T_in, win_size)
% MODULE2_LEE_FILTER_T
%   Polarimetric-preserving refined Lee filter applied to each element
%   of the coherency matrix independently (but using a single edge map
%   to preserve inter-channel phase relationships).
%
%   Reference: Lee et al., IEEE TGRS, 2008 (simplified implementation)

T = T_in;
fields = {'T11','T22','T33','T12','T13','T23'};

% Use Span as the reference image for edge/direction detection
Span = T_in.T11 + T_in.T22 + T_in.T33;
Span_dB = 10*log10(Span + eps);

% Estimate variance map from Span (log domain)
Span_mean = boxcar_avg(Span_dB, win_size);
Span_var  = boxcar_avg(Span_dB.^2, win_size) - Span_mean.^2;
Span_var  = max(Span_var, 0);

% Estimate noise variance (assumed ~ENL proportional)
% ENL estimate from uniform regions (use percentile)
noise_var = median(Span_var(:)) * 0.5;

% Weight function w (LMMSE weight)
W = Span_var ./ (Span_var + noise_var);
W = min(max(W, 0), 1);

% Apply to each T element
for fi = 1:length(fields)
    f = fields{fi};
    T_elem  = T_in.(f);
    T_local = boxcar_avg(T_elem, win_size);
    % LMMSE: filtered = mean + W*(original - mean)
    T.(f) = T_local + W .* (T_elem - T_local);
end

end


% =========================================================================
function out = boxcar_avg(x, win)
% Uniform boxcar spatial average using 2D convolution
    kernel = ones(win, win) / (win^2);
    if isreal(x)
        out = conv2(x, kernel, 'same');
    else
        out = conv2(real(x), kernel, 'same') + ...
           1j*conv2(imag(x), kernel, 'same');
    end
end
