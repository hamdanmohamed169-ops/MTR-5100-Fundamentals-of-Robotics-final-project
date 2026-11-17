function PL = fp_arm_plan(qs, qg, P, opt)
%FP_ARM_PLAN  Collision-free transport of the RCM unit by the 4-DOF arm
%   (Modules 8 and 9). No hand-picked via points:
%   1. C-space: 4 arm joints within limits; a configuration is free when the
%      arm + carried unit keep >= opt.dSafe mm from the patient (face SDF),
%      the body, the table and the camera stand (FP_ARM_CLEARANCE) AND the
%      arm does not hit itself: unit vs column / links, link vs link,
%      >= opt.dSelf mm (FP_ARM_SELFCLEAR).
%   2. RRT-Connect between the start and goal joints (FP_RRT_CONNECT).
%   3. Random shortcuts (FP_SHORTCUT), then corner rounding by Gaussian
%      smoothing of the densified path (the window is halved until the
%      smoothed path is collision-free again).
%   4. Time law with speed-and-separation monitoring: along the path the
%      speed limit is the largest v whose protective-stop distance
%      v*Tr + v^2/(2a) + margin fits in the local clearance (M9), capped by
%      P.arm.vmax (flange and RCM point) and P.arm.wmax (joints). The path
%      follows a smooth trapezoidal profile on this "time length".
%   5. Verification of the final trajectory with every mesh vertex.
%   PL.Q (4xN) joints at PL.t, PL.v / PL.vR flange / RCM-point speed (m/s),
%   PL.clr sampled-model clearance per frame (mm), PL.clrFull full-mesh
%   clearance (mm), PL.vlim speed limit, PL.ssm stop-distance margin (mm),
%   PL.raw / PL.short RRT / shortcut waypoints, PL.pF flange path (3xN, m), PL.stats.
if nargin < 4, opt = struct(); end
def = struct('dSafe', 12, 'dSelf', 10, 'seed', 7, 'dt', 0.04, 'step', 0.15, 'edgeRes', 0.01, 'maxIter', 5000, ...
    'nShort', 80, 'nShort2', 60, 'w', [1; 1; 1; 0.6], 'ramp', 0.8, 'verifyEvery', 2, 'tries', 3);
fn = fieldnames(def);
for i = 1:numel(fn), if ~isfield(opt, fn{i}), opt.(fn{i}) = def.(fn{i}); end, end
tic0 = tic;
CM = fp_arm_collision_model(P); ENV = fp_arm_env(P);
lo = P.arm.qmin; hi = P.arm.qmax;
clr = @(q) fp_arm_clearance(q, P, CM, ENV);
slf = @(q) fp_arm_selfclear(q, P, CM);
isFree = @(q) all(q >= lo & q <= hi) && slf(q) >= opt.dSelf && clr(q) >= opt.dSafe;
% adaptive edge check: no robot point moves more than sum_i |dq_i| r_i (r_i = reach
% beyond joint i), so a configuration with clearance c lets the edge advance by
% (c - dSafe)/L before the next check (never less than the resolution step)
rU = max(vecnorm(CM.unitFull)) + 60; r4 = P.arm.DH(4,3)*1e3 + rU; r3 = P.arm.DH(3,3)*1e3 + r4;
r2 = P.arm.DH(2,3)*1e3 + r3; rJ = [r2; r2; r3; r4];
opt.edgeFree = @(a, b) edgeAdaptive(a, b, clr, slf, lo, hi, opt.dSafe, opt.dSelf, rJ, opt.edgeRes);
rs = RandStream('mt19937ar', 'Seed', opt.seed); prev = RandStream.setGlobalStream(rs);
cln = onCleanup(@() RandStream.setGlobalStream(prev));
% best of opt.tries RRT-Connect runs (shortest weighted joint path after shortcutting);
% the second shortcut pass works on the densified path (finer shortcuts)
plen = @(Q) sum(sqrt(sum((opt.w(:).*diff(Q, 1, 2)).^2, 1)));
best = inf; tt = tic;
for tr = 1:opt.tries
    [raw_, I_] = fp_rrt_connect(isFree, qs(:), qg(:), lo, hi, opt);
    if ~I_.success, continue; end
    sh = shortcutA(raw_, opt.nShort, opt.edgeFree);
    if plen(sh) < best, best = plen(sh); raw = raw_; I = I_; short = sh; end
end
assert(isfinite(best), 'fp_arm_plan: RRT-Connect found no path in %d tries', opt.tries);
short = shortcutA(densify(short, 0.25), opt.nShort2, opt.edgeFree);   % finer shortcuts
tPlan = toc(tt); tt = tic;
% ---- densify + corner rounding ----
D = densify(short, 0.01);
win = 61; Qd = D; nOK = false;
while win >= 3
    Qd = smoothPath(D, win);
    nOK = true;
    for k = 2:size(Qd,2)-1, if ~isFree(Qd(:,k)), nOK = false; break; end, end
    if nOK, break; end
    win = 2*floor(win/4) + 1;
end
if ~nOK, Qd = D; end
% ---- clearance-dependent speed limit (speed and separation monitoring) ----
n = size(Qd,2); c = zeros(1,n); pF = zeros(3,n); pR = zeros(3,n);
for k = 1:n
    c(k) = clr(Qd(:,k));
    TF = fp_arm_fk(Qd(:,k), P); TR = TF*P.arm.T_40;
    pF(:,k) = TF(1:3,4); pR(:,k) = TR(1:3,4);
end
dLow = (c - CM.sampErr)/1e3;                                     % lower bound of the true clearance (m)
a = P.safe.arm_abrake; Tr = P.safe.arm_Tr; mg = P.safe.arm_margin;
vS = a*(-Tr + sqrt(Tr^2 + 2*max(dLow - mg, 0)/a));
vlim = max(min(P.arm.vmax, vS), 0.005);
vl = min(vlim(1:end-1), vlim(2:end));
dtau = max([vecnorm(diff(pF, 1, 2))./vl; vecnorm(diff(pR, 1, 2))./vl; max(abs(diff(Qd, 1, 2)), [], 1)/P.arm.wmax], [], 1);
tau = [0 cumsum(dtau)]; tauT = tau(end);
% smooth trapezoid on tau: d(tau)/dt ramps 0 -> 1 -> 0 with quintic smooth steps of length 'ramp'
Trp = min(opt.ramp, tauT); T = tauT + Trp;
t = 0:opt.dt:T; if t(end) < T - 1e-9, t = [t T]; end
tf = linspace(0, T, 4000);
sv = @(x) 10*x.^3 - 15*x.^4 + 6*x.^5;
vel = min(sv(min(tf/Trp, 1)), sv(min((T - tf)/Trp, 1)));
tauf = [0 cumsum((vel(1:end-1) + vel(2:end))/2.*diff(tf))]; tauf = tauf*tauT/tauf(end);
tauT_t = interp1(tf, tauf, t);
idx = interp1(tau, 1:n, min(tauT_t, tauT));
Q = interp1((1:n).', Qd.', idx).';
Q(:,end) = qg(:);
% ---- per-frame quantities + full-mesh verification ----
N = numel(t); v = zeros(1,N); vR = v; cl = v; vlimT = v; pFt = zeros(3,N); pRt = pFt;
for k = 1:N
    TF = fp_arm_fk(Q(:,k), P); TR = TF*P.arm.T_40; pFt(:,k) = TF(1:3,4); pRt(:,k) = TR(1:3,4);
end
v(2:end) = vecnorm(diff(pFt, 1, 2))/opt.dt; vR(2:end) = vecnorm(diff(pRt, 1, 2))/opt.dt;
cl = interp1(1:n, c, idx); vlimT = interp1(1:n, vlim, idx);   % sampled-model clearance / limit per frame
tSmooth = toc(tt); tt = tic;
kv = unique([1:opt.verifyEvery:N N]); cf = nan(1,N); body = zeros(1,N); obst = body;
cs = nan(1,N); spair = cell(1,N);
for k = kv
    [cf(k), inf_] = fp_arm_clearance(Q(:,k), P, CM, ENV, true); body(k) = inf_.body; obst(k) = inf_.obst;
    [cs(k), sinf] = fp_arm_selfclear(Q(:,k), P, CM, true); spair{k} = sinf.pair;
end
[csmin, ksmin] = min(cs);
tVerify = toc(tt);
vmaxT = max(v, vR);
sstop = (vmaxT*Tr + vmaxT.^2/(2*a))*1e3 + mg*1e3;                % protective-stop distance + margin (mm)
ssm = cf - sstop;                                                % >= 0 -> the arm can always stop before contact
[cmin, kmin] = min(cf);
PL = struct('Q', Q, 't', t, 'T', T, 'v', v, 'vR', vR, 'clr', cl, 'clrFull', cf, 'vlim', vlimT, ...
    'ssm', ssm, 'clrSelf', cs, 'raw', raw, 'short', short, 'dense', Qd, 'pF', pFt, 'pR', pRt, 'bodyMin', body, 'obstMin', obst);
PL.stats = struct('iter', I.iter, 'nodes', I.nodes, 'checks', I.checks, 'nRaw', size(raw,2), 'nShort', size(short,2), ...
    'smoothWin', win*nOK, 'minClr', cmin, 'minBody', CM.names{body(kmin)}, 'minObst', ENV.names{obst(kmin)}, ...
    'minSSM', min(ssm), 'pathLenF', sum(vecnorm(diff(pFt, 1, 2))), 'time', toc(tic0), 'dSafe', opt.dSafe, ...
    'sampErr', CM.sampErr, 'okClr', cmin > 0, 'okSSM', min(ssm) >= 0, 'tPlan', tPlan, 'tSmooth', tSmooth, 'tVerify', tVerify, 'minSelf', csmin, 'minSelfPair', spair{ksmin}, ...
    'dSelf', opt.dSelf, 'okSelf', csmin > 0);
end

function [ok, nc] = edgeAdaptive(a, b, clr, slf, lo, hi, dSafe, dSelf, rJ, res)
% collision check of the straight joint-space edge a -> b with clearance-based steps
% (self distances can shrink at most twice as fast: both bodies move)
dq = b - a; L = sum(abs(dq).*rJ); dsMin = res/max(max(abs(dq)), eps);
s = 0; nc = 0; ok = true;
while true
    q = a + s*dq;
    if any(q < lo | q > hi), ok = false; return; end
    cS = slf(q); nc = nc + 1;
    if cS < dSelf, ok = false; return; end
    c = clr(q);
    if c < dSafe, ok = false; return; end
    if s >= 1, return; end
    s = min(1, s + max(min((c - dSafe)/max(L, eps), (cS - dSelf)/max(2*L, eps)), dsMin));
end
end

function path = shortcutA(path, nTry, edgeFree)
% random shortcuts (FP_SHORTCUT) with the adaptive edge check
for k = 1:nTry
    K = size(path,2);
    if K <= 2, break; end
    ij = sort(randperm(K, 2)); i = ij(1); j = ij(2);
    if j - i < 2, continue; end
    if edgeFree(path(:,i), path(:,j)), path = [path(:,1:i) path(:,j:end)]; end
end
end

function D = densify(Q, ds)
% insert points so that consecutive waypoints differ by <= ds (rad, max norm)
D = Q(:,1);
for k = 1:size(Q,2)-1
    m = max(1, ceil(max(abs(Q(:,k+1) - Q(:,k)))/ds));
    D = [D Q(:,k) + (Q(:,k+1) - Q(:,k))*((1:m)/m)]; %#ok<AGROW>
end
end

function S = smoothPath(D, win)
% Gaussian smoothing of a dense path with fixed end points (end values padded)
h = ceil(win/2); g = exp(-0.5*((-h:h)/(win/5)).^2); g = g/sum(g);
Pd = [repmat(D(:,1), 1, h) D repmat(D(:,end), 1, h)];
S = zeros(size(D));
for i = 1:size(D,1), S(i,:) = conv(Pd(i,:), g, 'valid'); end
S(:,1) = D(:,1); S(:,end) = D(:,end);
end
