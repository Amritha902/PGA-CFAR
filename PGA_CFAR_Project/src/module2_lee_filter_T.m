function T = module2_lee_filter_T(T_in, win_size)
% MODULE2_LEE_FILTER_T — standalone wrapper (called from main)
% Delegates to the implementation inside module2_build_coherency.m
% This file exists only so MATLAB can resolve the function name directly.

fields = {'T11','T22','T33','T12','T13','T23'};
T = T_in;

Span = T_in.T11 + T_in.T22 + T_in.T33;
Span_dB = 10*log10(Span + eps);

kernel = ones(win_size, win_size) / (win_size^2);
Span_mean = conv2(Span_dB, kernel, 'same');
Span_var  = conv2(Span_dB.^2, kernel, 'same') - Span_mean.^2;
Span_var  = max(Span_var, 0);

noise_var = median(Span_var(:)) * 0.5;
W = Span_var ./ (Span_var + noise_var + eps);
W = min(max(W, 0), 1);

for fi = 1:length(fields)
    f = fields{fi};
    T_elem = T_in.(f);
    if isreal(T_elem)
        T_local = conv2(T_elem, kernel, 'same');
    else
        T_local = conv2(real(T_elem), kernel, 'same') + ...
               1j*conv2(imag(T_elem), kernel, 'same');
    end
    T.(f) = T_local + W .* (T_elem - T_local);
end
end
