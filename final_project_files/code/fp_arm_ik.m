function [q, info] = fp_arm_ik(X, P, qprev)
%FP_ARM_IK  Closed-form inverse kinematics of the 4-DOF positioning arm
%   (base yaw + 3-link planar arm, Module 5).
%   X     : desired flange pose, either a 4x4 transform in {E} (metres), or
%           a task vector [x; y; z; phi] = flange position in {E} + pitch
%           angle phi of the flange x-axis above the horizontal
%   qprev : previous joint vector -> the closest valid branch is returned
%   info.Q      : the 4 candidate solutions (4x4, columns)
%   info.ok     : candidate reaches the pose and respects the joint limits
%   info.err    : [position (m) orientation (rad)] residual of FK(q) vs X
%   info.branch : {'front'|'back', 'elbow up'|'elbow down'}
%
%   Step 1 - base yaw:   q1 = atan2(y_b, x_b)   (or q1 + pi, "back")
%   Step 2 - wrist point in the arm plane:  (r_w, s_w) = (r, s) - a4 (cos phi, sin phi)
%   Step 3 - planar 2-link IK (law of cosines):
%            cos q3 = (r_w^2 + s_w^2 - a2^2 - a3^2) / (2 a2 a3),  q3 = +-acos(.)
%            q2 = atan2(s_w, r_w) - atan2(a3 sin q3, a2 + a3 cos q3)
%   Step 4 - wrist:      q4 = phi - q2 - q3
if nargin < 3, qprev = zeros(4,1); end
A = P.arm;
if isequal(size(X), [4 4])
    Tg = X; pE = X(1:3,4); xE = X(1:3,1);
else
    pE = X(1:3); Tg = [];
end
Tb = inv(A.T_Eb);
pb = Tb(1:3,:)*[pE; 1];                       % flange position in the arm base frame
Q = nan(4,4); lab = cell(2,4); col = 0;
for sh = [1 -1]                               % front / back
    q1 = atan2(sh*pb(2), sh*pb(1));
    hdir = [cos(q1); sin(q1); 0];             % horizontal radial direction of the arm plane
    r = hdir.'*pb(1:3); s = pb(3) - A.d1;
    if isempty(Tg)
        phi = X(4);
    else
        xb = Tb(1:3,1:3)*xE;                  % flange x-axis in the base frame
        phi = atan2(xb(3), hdir.'*xb);
    end
    rw = r - A.a4*cos(phi); sw = s - A.a4*sin(phi);
    D = (rw^2 + sw^2 - A.a2^2 - A.a3^2)/(2*A.a2*A.a3);
    for el = [-1 1]                           % elbow up (q3 < 0) / elbow down
        col = col + 1;
        lab(:,col) = {pick(sh > 0, 'front', 'back'); pick(el < 0, 'elbow up', 'elbow down')};
        if abs(D) > 1, continue; end          % wrist point out of reach
        q3 = el*acos(D);
        q2 = atan2(sw, rw) - atan2(A.a3*sin(q3), A.a2 + A.a3*cos(q3));
        q4 = phi - q2 - q3;
        Q(:,col) = [q1; q2; q3; q4];
    end
end
ok = false(1,4); dist = inf(1,4); err = inf(4,2);
for c = 1:4
    if any(isnan(Q(:,c))), continue; end
    Q(:,c) = qprev + atan2(sin(Q(:,c) - qprev), cos(Q(:,c) - qprev));   % unwrap towards qprev
    Tc = fp_arm_fk(Q(:,c), P);
    err(c,1) = norm(Tc(1:3,4) - pE);
    if ~isempty(Tg)
        dR = Tc(1:3,1:3).'*Tg(1:3,1:3);
        err(c,2) = acos(min(1, max(-1, (trace(dR) - 1)/2)));
    else
        err(c,2) = 0;
    end
    ok(c) = all(Q(:,c) >= A.qmin - 1e-9 & Q(:,c) <= A.qmax + 1e-9) && err(c,1) < 1e-9 && err(c,2) < 1e-6;
    if ok(c), dist(c) = norm(Q(:,c) - qprev); end
end
[~, best] = min(dist);
if ~ok(best)
    q = qprev; info = struct('Q', Q, 'ok', ok, 'best', 0, 'branch', {{'-', '-'}}, 'err', [inf inf]);
    return
end
q = Q(:,best);
info = struct('Q', Q, 'ok', ok, 'best', best, 'branch', {lab(:,best).'}, 'err', err(best,:));
end

function s = pick(c, a, b)
if c, s = a; else, s = b; end
end
