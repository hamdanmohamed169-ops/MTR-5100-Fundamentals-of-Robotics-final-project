function [p, u] = fp_fk_geom(q)
%FP_FK_GEOM  Direct geometric (spherical-coordinate) FK, independent of DH.
%   The needle passes through the RCM (origin of {0}); its unit direction is
%   u = [cos q1 sin q2; sin q1 sin q2; cos q2] and the tip is p = q3*u.
u = [cos(q(1))*sin(q(2)); sin(q(1))*sin(q(2)); cos(q(2))];
p = q(3)*u;
end
