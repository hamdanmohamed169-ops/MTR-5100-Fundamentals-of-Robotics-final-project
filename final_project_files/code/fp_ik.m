function [q, info] = fp_ik(p, P, qprev)
%FP_IK  Analytical inverse kinematics of the 3-DOF RCM robot (Module 5).
%   p     : desired needle-tip position in {0} (3x1, m)
%   qprev : previous configuration used to pick the branch (optional)
%
%   Spherical coordinates of p (R-R-P robot):
%       q3 = +-|p|,   q2 = acos(u_z),   q1 = atan2(u_y, u_x),   u = p/q3
%   Four algebraic candidates exist (like elbow-up/down of the 2R arm):
%     A: q2 in (0,pi),  q3 = +|p|       B: q1+pi, -q2,  q3 = +|p|
%     C: needle reversed, q3 = -|p|     D: C with q1+pi, -q2
%   Candidates are filtered by joint limits and the one closest to qprev
%   (wrapped angle distance) is returned. info lists all candidates.
if nargin < 3 || isempty(qprev), qprev = P.q_home; end
d = norm(p);
info.singular = d < 1e-9;
cand = zeros(3,4);
if info.singular
    % tip at the RCM: q1,q2 undetermined -> keep previous orientation
    cand(:,1) = [qprev(1); qprev(2); 0];
    ncand = 1;
else
    k = 0;
    for s3 = [1 -1]
        u  = s3*p/d;
        q2 = acos(max(-1,min(1,u(3))));
        q1 = atan2(u(2), u(1));
        k = k + 1; cand(:,k) = [q1;  q2; s3*d];
        k = k + 1; cand(:,k) = [wrapToPiLocal(q1+pi); -q2; s3*d];
    end
    ncand = 4;
end
cand = cand(:,1:ncand);
ok = all(cand >= P.qmin - 1e-12 & cand <= P.qmax + 1e-12, 1);
info.cand = cand;  info.within_limits = ok;
% FK reconstruction error of every candidate
info.err = zeros(1,ncand);
for k = 1:ncand
    info.err(k) = norm(fp_fk_geom(cand(:,k)) - p);
end
if ~any(ok)
    q = []; info.status = 'no solution within joint limits';
    return
end
idx = find(ok);
dist = zeros(size(idx));
for k = 1:numel(idx)
    dq = cand(:,idx(k)) - qprev;
    dq(1:2) = atan2(sin(dq(1:2)), cos(dq(1:2)));
    dist(k) = norm(dq(1:2)) + 100*abs(dq(3));   % 100 rad/m weighting
end
[~, j] = min(dist);
q = cand(:,idx(j));
info.status = 'ok';
info.selected = idx(j);
end

function a = wrapToPiLocal(a)
a = atan2(sin(a), cos(a));
end
