function [dmin, info] = fp_arm_selfclear(qa, P, CM, full)
%FP_ARM_SELFCLEAR  Self-collision clearance (mm) of the 4-DOF arm and the RCM
%   unit it carries (robot at q_home), for arm joints qa (4x1, rad).
%   The arm bodies are capsules sized from the drawn geometry (FP_SCENE /
%   buildArm): column C0, shoulder housing C1h, upper arm C2, elbow housing
%   C2h, forearm C3, wrist housing C3h, wrist link + flange C4. The unit is
%   the point cloud of FP_ARM_COLLISION_MODEL (sampled, or every vertex when
%   full = true). Only NON-adjacent pairs are checked; segments are trimmed
%   near a shared joint so that the joint housings themselves do not count.
%   info.pair : name of the closest pair
if nargin < 4, full = false; end
[TE4, Ta] = fp_arm_fk(qa, P);
mm = 1e3; Tb = P.arm.T_Eb;
B = Tb(1:3,4)*mm; zb = Tb(1:3,3);
S = Ta(1:3,4,1)*mm; E = Ta(1:3,4,2)*mm; W = Ta(1:3,4,3)*mm; F = TE4(1:3,4)*mm;
z1 = Ta(1:3,3,1); z2 = Ta(1:3,3,2); z3 = Ta(1:3,3,3);
u2 = (E - S)/norm(E - S); u3 = (W - E)/norm(W - E);
d1 = P.arm.d1*mm;
C = struct('n', {}, 'a', {}, 'b', {}, 'r', {});
C(1) = cap('column',          B - 90*zb,  B + (d1 - 80)*zb, 44);
C(2) = cap('shoulder housing', S - 56*z1, S + 56*z1, 42);
C(3) = cap('upper arm',        S + 90*u2, E, 30);
C(4) = cap('elbow housing',    E - 44*z2, E + 44*z2, 33);
C(5) = cap('forearm',          E + 75*u3, W, 25);
C(6) = cap('wrist housing',    W - 34*z3, W + 34*z3, 26);
C(7) = cap('wrist link',       W, F, 24);
C(8) = cap('upper arm (elbow side)', S, E - 75*u2, 30);
pairs = [1 3; 1 4; 1 5; 1 6; 1 7; 2 5; 2 6; 2 7; 8 5; 3 6; 3 7; 4 7];
dmin = inf; info = struct('pair', '');
for k = 1:size(pairs,1)
    a = C(pairs(k,1)); b = C(pairs(k,2));
    d = segSeg(a.a, a.b, b.a, b.b) - a.r - b.r;
    if d < dmin, dmin = d; info.pair = [a.n ' / ' b.n]; end
end
% RCM unit points against every arm body except the wrist (it carries the unit)
if full, U = CM.unitFull; else, U = CM.unit; end
Uw = TE4(1:3,1:3)*U + TE4(1:3,4)*mm;
for k = [1 2 3 4 5]
    c = C(k); ab = c.b - c.a; t = min(max((ab.'*(Uw - c.a))/(ab.'*ab), 0), 1);
    d = min(sqrt(sum((Uw - c.a - ab*t).^2, 1))) - c.r;
    if d < dmin, dmin = d; info.pair = ['RCM unit / ' c.n]; end
end
end

function c = cap(n, a, b, r)
c = struct('n', n, 'a', a, 'b', b, 'r', r);
end

function d = segSeg(p1, q1, p2, q2)
% shortest distance between segments p1-q1 and p2-q2
d1 = q1 - p1; d2 = q2 - p2; r = p1 - p2;
a = d1.'*d1; e = d2.'*d2; f = d2.'*r; c = d1.'*r; b = d1.'*d2;
den = a*e - b^2;
if den > 1e-12, s = min(max((b*f - c*e)/den, 0), 1); else, s = 0; end
t = (b*s + f)/e;
if t < 0, t = 0; s = min(max(-c/a, 0), 1);
elseif t > 1, t = 1; s = min(max((b - c)/a, 0), 1); end
d = norm(p1 + d1*s - (p2 + d2*t));
end
