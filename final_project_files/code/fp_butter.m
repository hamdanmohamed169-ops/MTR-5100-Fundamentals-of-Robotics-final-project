function [sos, H, f] = fp_butter(order, fc, fs, fgrid)
%FP_BUTTER  Digital Butterworth low-pass filter designed by hand with the
%   bilinear transform (no Signal Processing Toolbox needed).
%   sos : rows [b0 b1 b2 1 a1 a2] of cascaded 2nd-order sections
%         (order 1 uses b2 = a2 = 0).
%   H   : complex frequency response evaluated at frequencies fgrid (Hz).
%   Pre-warping: K = tan(pi fc / fs). Each 2nd-order section has quality
%   factor Q_k = 1/(2 cos((2k-1) pi/(2n))).
K = tan(pi*fc/fs);
sos = [];
if mod(order,2) == 1
    sos = [K/(1+K) K/(1+K) 0 1 (K-1)/(K+1) 0];
end
for k = 1:floor(order/2)
    Q = 1/(2*cos((2*k-1)*pi/(2*order)));
    nrm = 1 + K/Q + K^2;
    sos = [sos; K^2/nrm 2*K^2/nrm K^2/nrm 1 2*(K^2-1)/nrm (1-K/Q+K^2)/nrm]; %#ok<AGROW>
end
if nargout > 1
    if nargin < 4, fgrid = logspace(-1, log10(fs/2*0.99), 400); end
    z = exp(1i*2*pi*fgrid/fs);
    H = ones(size(z));
    for r = 1:size(sos,1)
        H = H .* (sos(r,1) + sos(r,2)*z.^-1 + sos(r,3)*z.^-2) ./ ...
                 (1 + sos(r,5)*z.^-1 + sos(r,6)*z.^-2);
    end
    f = fgrid;
end
end
