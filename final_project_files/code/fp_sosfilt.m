function y = fp_sosfilt(sos, x)
%FP_SOSFILT  Apply cascaded sections from FP_BUTTER using base-MATLAB filter().
%   Initial conditions are set to the steady state of the first sample so
%   the filter does not produce a start-up transient (signal starts at rest).
y = x(:);
for r = 1:size(sos,1)
    b = sos(r,1:3); a = sos(r,4:6);
    % steady-state initial condition for constant input y(1)
    zi = steady_ic(b, a) * y(1);
    y = filter(b, a, y, zi);
end
y = reshape(y, size(x));
end

function zi = steady_ic(b, a)
% Direct-form-II-transposed states for unit DC input with unit DC gain
% (filter() uses 2 states for a 3-coefficient section).
g  = sum(b)/sum(a);           % DC gain (=1 for low-pass)
z2 = b(3) - a(3)*g;
z1 = b(2) + z2 - a(2)*g;
zi = [z1; z2];
end
