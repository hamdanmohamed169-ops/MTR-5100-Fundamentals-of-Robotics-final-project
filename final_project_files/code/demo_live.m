function demo_live(varargin)
%DEMO_LIVE  Live 3-D prototype demonstration of the complete intravitreal
%   injection procedure:
%     - a stereo camera locates the eye (triangulation + circle fit),
%     - a 4-DOF positioning arm carries the injection unit to the eye along a
%       collision-free RRT-Connect path (patient face SDF, body, table and
%       camera stand as obstacles, speed limited by the stopping distance)
%       and docks it on the camera estimate with closed-form inverse kinematics,
%     - the 3-DOF RCM robot (yaw, pitch, insertion) performs the injection,
%     - the arm returns to its exact starting pose with inverse kinematics.
%
%   demo_live              - play the animation (about 1.5 minutes)
%   demo_live('record')    - play and also save ../results/demo_procedure.mp4
%   demo_live('fast')      - play without real-time pacing
%   demo_live('snapshots') - also save one PNG per step to ../figures (D01..D15)
%
%   Everything shown is computed with the project models (no pre-recorded
%   data): camera registration (fp_eye_register), arm FK/IK (fp_arm_fk /
%   fp_arm_ik) along Cartesian paths, RCM DH kinematics + IK (fp_fk/fp_ik),
%   trajectories (fp_traj), needle-force model (fp_needle_force),
%   A*-planned re-orientation (fp_astar + fp_shortcut), IOP-feedback
%   injection (fp_iop_sim) and the eye-twitch safety layer (fp_twitch_sim).
record = any(strcmpi(varargin, 'record'));
fast   = any(strcmpi(varargin, 'fast'));
snaps  = any(strcmpi(varargin, 'snapshots'));
dry    = any(strcmpi(varargin, 'dryrun'));        % build the timeline only (testing)
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
P = fp_params(); rng(P.seed);
mm = 1e3; um = 1e6; dt = 0.04;

%% ---------------- 1. Build the procedure timeline ----------------------
fprintf('Preparing the procedure (camera, arm IK, RCM kinematics, planning, fluidics)...\n');
phaseName = {'1  STEREO CAMERA LOCATES THE EYE', ...
    '2  LIMBUS CIRCLE FIT  ->  EYE CENTRE AND ENTRY POINT', ...
    '3  INVERSE KINEMATICS OF THE DOCKING POSE (4-DOF ARM)', ...
    '4  TRANSPORT ALONG THE RRT-CONNECT PATH', ...
    '5  DOCKING:  RCM POINT ONTO THE MEASURED ENTRY POINT', ...
    '6  ALIGN NEEDLE WITH THE PLANNED ENTRY AXIS', '7  APPROACH THE SCLERA (PARS PLANA)', ...
    '8  INSERTION THROUGH THE SCLERA', '9  DRUG INJECTION AT SITE A (IOP FEEDBACK)', ...
    '10  RE-ORIENTATION A -> B ABOUT THE RCM (A* PATH)', '11  DRUG INJECTION AT SITE B (IOP FEEDBACK)', ...
    '12  SAFETY TEST: SUDDEN EYE TWITCH', '13  RETRACT THE NEEDLE', ...
    '14  UNDOCK ALONG THE INSERTION AXIS', '15  RETURN TO THE START POSE (RRT-CONNECT)'};
nPh = numel(phaseName);
TL = struct('q',{},'qa',{},'ph',{},'msg',{},'F',{},'iop',{},'V',{},'th',{},'tp',{},'pace',{},'arm',{},'cam',{});
tproc = 0;

% ---- camera registration ----
REG = fp_eye_register(P, P.seed);
nP = P.cam.nPts;

% ---- arm key poses (task space: flange position + pitch angle phi) ----
x0 = P.T_E0(1:3,1);
phiDock = taskOf(P.T_E0/P.arm.T_40, P); phiDock = phiDock(4);   % pitch of the flange when docked
xDock = taskForRCM(REG.entry, phiDock, P);          % flange task pose that puts the RCM on the MEASURED entry
xPark = [0.34; -0.36; 0.05; -pi/2];                  % start: parked low beside the patient's chest (surgeon's view free)
[qPark, Ipark0] = fp_arm_ik(xPark, P, deg2rad([60; 80; -100; -70]));
[qDockSol, Idock] = fp_arm_ik(xDock, P, qPark);
TE0d = fp_arm_fk(qDockSol, P)*P.arm.T_40;           % docked RCM frame (insertion axis = its x-axis)
xPre  = xDock; xPre(1:3) = xPre(1:3) - 0.040*TE0d(1:3,1);   % 40 mm back along the insertion axis
qPre  = fp_arm_ik(xPre, P, qDockSol);
fprintf('Planning the arm transport (RRT-Connect in the 4-D joint space, patient SDF as obstacle)...\n');
planOpt = struct('dSafe', 25, 'dt', dt, 'verifyEvery', 4);
PLo = fp_arm_plan(qPark, qPre, P, planOpt);
assert(Ipark0.best > 0 && Idock.best > 0, 'arm pose not reachable');
none = armInfo(qPark, xPark, P, 0, Ipark0); none.br = 'parked';
camOff = struct('n', 0, 'fit', false, 'rays', false);

% --- phases 1-2: camera registration (robot parked) ---
for k = 1:round(3.2/dt)
    c = camOff; c.n = min(nP, ceil(nP*k/round(2.6/dt))); c.rays = true;
    TL(end+1) = frame(P.q_home, qPark, 1, sprintf(['Both cameras detect limbus edge points (%d / %d)  ->  ' ...
        'linear triangulation (DLT, SVD) gives their 3-D positions'], c.n, nP), 0, P.iop.P0, 0, 0, 0, dt, none, c); %#ok<AGROW>
end
for k = 1:round(3.2/dt)
    c = camOff; c.n = nP; c.fit = true;
    TL(end+1) = frame(P.q_home, qPark, 2, sprintf(['Plane + circle fit: limbus r = %.2f mm  ->  eye centre error %.0f um,  ' ...
        'entry-point error %.0f um (limit %.0f um)'], REG.limbus_r*mm, REG.err_centre*um, REG.err_entry*um, P.cam.req*um), ...
        0, P.iop.P0, 0, 0, 0, dt, none, c); %#ok<AGROW>
end
camKeep = camOff; camKeep.fit = true;               % keep the measured entry point on screen
% --- phase 3: IK of the docking pose (ghost preview) ---
for k = 1:round(2.4/dt)
    a = armInfo(qDockSol, xDock, P, 0, Idock);       % residual / |det J| of the IK solution being previewed
    TL(end+1) = frame(P.q_home, qPark, 3, sprintf(['Closed-form IK: q1 = atan2, planar 2-link law of cosines, q4 = phi - q2 - q3  ->  ' ...
        '%d of 4 solutions valid, chosen %s / %s'], nnz(Idock.ok), Idock.branch{:}), 0, P.iop.P0, 0, 0, tproc, dt, a, camKeep); %#ok<AGROW>
end
% --- phase 4: transport along the RRT-Connect path (no hand-picked via points) ---
[QA, AI, T4] = planFrames(PLo, xPre, P);
st4 = PLo.stats;
for k = 1:size(QA,2)
    TL(end+1) = frame(P.q_home, QA(:,k), 4, sprintf(['RRT-Connect: %d iter, %d -> %d waypoints  |  speed %5.1f mm/s ' ...
        '(SSM limit %5.1f)  |  clearance %5.1f mm (min %4.1f)'], st4.iter, st4.nRaw, st4.nShort, AI(k).v*mm, ...
        PLo.vlim(k)*mm, PLo.clr(k), st4.minClr), 0, P.iop.P0, 0, 0, tproc + (k-1)*dt, dt, AI(k), camKeep); %#ok<AGROW>
end
tproc = tproc + T4; pathOut = AI(1).path; rawOut = rawTrace(PLo, P);
% --- phase 5: docking, straight line along the insertion axis ---
[QA, AI, T5] = armPath(xPre, xDock, [], QA(:,end), P, P.arm.vdock, dt);
for k = 1:size(QA,2)
    TL(end+1) = frame(P.q_home, QA(:,k), 5, sprintf('Straight-line docking at <= %g mm/s   RCM-to-measured-entry distance %5.2f mm', ...
        P.arm.vdock*mm, AI(k).d*mm), 0, P.iop.P0, 0, 0, tproc + (k-1)*dt, dt, AI(k), camKeep); %#ok<AGROW>
end
tproc = tproc + T5; qDock = QA(:,end); pathDock = AI(1).path;
TE0dock = fp_arm_fk(qDock, P)*P.arm.T_40;
dockErrEst = norm(TE0dock(1:3,4) - REG.entry); dockErrTrue = norm(TE0dock(1:3,4) - P.eye.pE);
lockedArm = armInfo(qDock, xDock, P, 0, Idock); lockedArm.br = 'locked';
for k = 1:round(1.6/dt)
    TL(end+1) = frame(P.q_home, qDock, 5, sprintf(['Docked: RCM on the measured entry (%.0e mm); true entry point %.0f um away ' ...
        '(camera error)  ->  arm brakes locked'], dockErrEst*mm, dockErrTrue*um), 0, P.iop.P0, 0, 0, tproc, dt, lockedArm, camOff); %#ok<AGROW>
end

% ---- RCM procedure (3-DOF robot, arm locked at qDock) ----
qA = [deg2rad(-20); deg2rad(100); 7e-3];
qB = [deg2rad( 20); deg2rad(100); 7e-3];
opt = struct('vmax', P.vmax, 'amax', P.amax);
RI = fp_iop_sim('feedback', P, struct('Tend', 600, 'margin', 1, 'Kp', 0.5, 'Qcap', 5, 'noise', false));
iopNow = P.iop.P0; Vnow = 0;
% phases 6, 7: trapezoidal joint moves
qAout = [qA(1:2); -10e-3]; qAsurf = [qA(1:2); -0.5e-3];
segs = {P.q_home, qAout, 6; qAout, qAsurf, 7};
for s = 1:2
    [~,~,~,T] = fp_traj('trapezoid', segs{s,1}, segs{s,2}, 0, opt);
    t = 0:dt:T; if t(end) < T, t = [t T]; end %#ok<AGROW>
    Q = fp_traj('trapezoid', segs{s,1}, segs{s,2}, t, opt);
    qik = fp_ik(fp_fk_geom(segs{s,2}), P, segs{s,1});
    for k = 1:numel(t)
        TL(end+1) = frame(Q(:,k), qDock, segs{s,3}, sprintf('RCM robot IK target  q = [%.1f deg, %.1f deg, %.1f mm]', ...
            rad2deg(qik(1:2)), qik(3)*mm), 0, iopNow, Vnow, 0, tproc + t(k), dt, lockedArm, camOff); %#ok<AGROW>
    end
    tproc = tproc + T;
end
% phase 8: quintic insertion with the needle-force model
Tin = 6; tf = 0:1e-3:Tin;
[Q, Qd] = fp_traj('quintic', qAsurf, qA, tf, struct('T', Tin));
st.punctured = false; msg = 'Quintic insertion (joint 3), needle-force model active'; Fk = zeros(size(tf)); punct = false(size(tf));
for k = 1:numel(tf)
    [Fk(k), st] = fp_needle_force(Q(3,k), Qd(3,k), st, P); punct(k) = st.punctured;
end
kp = find(punct, 1);
for k = 1:40:numel(tf)
    if k >= kp, msg = sprintf('Sclera punctured at q_3 = %.2f mm  ->  force drop detected (3 ms)', Q(3,kp)*mm); end
    TL(end+1) = frame(Q(:,k), qDock, 8, msg, Fk(k), iopNow, Vnow, 0, tproc + tf(k), dt, lockedArm, camOff); %#ok<AGROW>
end
tproc = tproc + Tin; Fhold = Fk(end);
% phase 9: injection A (first 25 uL), time-lapse
iA = find(RI.V >= P.iop.Vdose/2, 1);
idx = unique(round(linspace(1, iA, 45)));
for k = idx
    TL(end+1) = frame(qA, qDock, 9, sprintf('Injected %.1f / %.0f uL   IOP %.1f mmHg (limit %g)   time-lapse x%.0f', ...
        RI.V(k), P.iop.Vdose, RI.P(k), P.iop.Plim, RI.t(iA)/(45*0.06)), ...
        Fhold, RI.P(k), RI.V(k), 0, tproc + RI.t(k), 0.06, lockedArm, camOff); %#ok<AGROW>
end
tproc = tproc + RI.t(iA); iopNow = RI.P(iA); Vnow = RI.V(iA);
% phase 10: A* + shortcut planned re-orientation
d = qA(3); q1g = -40:1:40; q2g = 50:1:130;
occ = false(numel(q2g), numel(q1g));
for i = 1:numel(q2g)
    for j = 1:numel(q1g)
        occ(i,j) = ~fp_lens_check([deg2rad(q1g(j)); deg2rad(q2g(i)); d], P);
    end
end
isFree = @(qq) fp_lens_check([deg2rad(qq(1)); deg2rad(qq(2)); d], P);
pth = fp_astar(occ, [find(q2g==100) find(q1g==-20)], [find(q2g==100) find(q1g==20)]);
wp = fp_shortcut(isFree, [q1g(pth(:,2)); q2g(pth(:,1))], 300, 0.5);
sC = [0 cumsum(vecnorm(diff(wp,1,2)))]; Tp = 8.4;
t = 0:dt:Tp; s = sC(end)*(10*(t/Tp).^3 - 15*(t/Tp).^4 + 6*(t/Tp).^5);
q12 = deg2rad([pchip(sC, wp(1,:), s); pchip(sC, wp(2,:), s)]);
for k = 1:numel(t)
    TL(end+1) = frame([q12(:,k); d], qDock, 10, 'Needle pivots about the RCM: zero lateral motion at the sclera, no second puncture', ...
        Fhold, iopNow, Vnow, 0, tproc + t(k), 0.03, lockedArm, camOff); %#ok<AGROW>
end
tproc = tproc + Tp;
% phase 11: injection B (remaining dose)
iB = find(RI.V >= P.iop.Vdose - 1e-9, 1);
idx = unique(round(linspace(iA, iB, 50)));
for k = idx
    TL(end+1) = frame(qB, qDock, 11, sprintf('Injected %.1f / %.0f uL   IOP %.1f mmHg (limit %g)   time-lapse x%.0f', ...
        RI.V(k), P.iop.Vdose, RI.P(k), P.iop.Plim, (RI.t(iB)-RI.t(iA))/(50*0.06)), ...
        Fhold, RI.P(k), RI.V(k), 0, tproc + RI.t(k) - RI.t(iA), 0.06, lockedArm, camOff); %#ok<AGROW>
end
tproc = tproc + RI.t(iB) - RI.t(iA); Vnow = P.iop.Vdose; iopNow = RI.P(iB);
% phase 12: twitch with spring quick-release (slow motion x10)
RT = fp_twitch_sim(P.safe.twitch_T, qB, P, 'spring', P.seed + 31);
for k = 1:20:find(RT.t <= 0.30, 1, 'last')
    qk = [qB(1:2); RT.q4(k)];
    if RT.t(k) < RT.tdet
        msg = sprintf('Eye rotating: %.1f deg  (slow motion x10)', rad2deg(RT.th(k)));
    elseif isnan(RT.texit) || RT.t(k) < RT.texit
        msg = sprintf('TWITCH DETECTED by the camera at %.0f ms -> spring quick-release retraction', RT.tdet*1e3);
    else
        msg = sprintf('Needle out at %.0f ms, lateral scleral displacement %.2f mm (limit 0.5)', RT.texit*1e3, RT.lat_exit*mm);
    end
    F = P.needle.f_fric*max(qk(3),0);
    TL(end+1) = frame(qk, qDock, 12, msg, F, iopNow, Vnow, RT.th(k), tproc + RT.t(k), 0.02, lockedArm, camOff); %#ok<AGROW>
end
tproc = tproc + 0.30;
% phase 13: eye re-fixates, needle back to the RCM home configuration
i30 = find(RT.t <= 0.30, 1, 'last');
q8a = [qB(1:2); RT.q4(i30)]; q8b = [qB(1:2); -10e-3]; th0 = RT.th(i30);
[~,~,~,T] = fp_traj('trapezoid', q8a, q8b, 0, opt); T = max(T, 1);
t = 0:dt:T; Q = fp_traj('trapezoid', q8a, q8b, t, struct('vmax',P.vmax,'amax',P.amax,'T',T));
for k = 1:numel(t)
    TL(end+1) = frame(Q(:,k), qDock, 13, 'Patient re-fixates; needle withdrawn to the safe standby distance', 0, iopNow, Vnow, ...
        th0*(1 - min(t(k)/T,1)), tproc + t(k), dt, lockedArm, camOff); %#ok<AGROW>
end
tproc = tproc + T;
[~,~,~,T] = fp_traj('trapezoid', q8b, P.q_home, 0, opt);
t = 0:dt:T; if t(end) < T, t = [t T]; end
Q = fp_traj('trapezoid', q8b, P.q_home, t, opt);
for k = 1:numel(t)
    TL(end+1) = frame(Q(:,k), qDock, 13, 'Procedure complete: 50 uL delivered, peak IOP below 25 mmHg', 0, iopNow, Vnow, 0, ...
        tproc + t(k), dt, lockedArm, camOff); %#ok<AGROW>
end
tproc = tproc + t(end);

% ---- arm returns: undock along the axis, then IK back to the start pose ----
[QA, AI, T14] = armPath(xDock, xPre, [], qDock, P, P.arm.vdock, dt);
for k = 1:size(QA,2)
    TL(end+1) = frame(P.q_home, QA(:,k), 14, sprintf('Brakes released; straight-line retreat along the insertion axis   %5.1f mm from the eye', ...
        (0.040 - AI(k).d)*mm), 0, iopNow, Vnow, 0, tproc + (k-1)*dt, dt, AI(k), camOff); %#ok<AGROW>
end
tproc = tproc + T14;
[~, Ipark] = fp_arm_ik(xPark, P, QA(:,end));
for k = 1:round(1.6/dt)
    a = armInfo(QA(:,end), xPark, P, 0, Ipark);
    TL(end+1) = frame(P.q_home, QA(:,end), 15, sprintf('IK of the start pose: %d of 4 solutions valid, chosen %s / %s (closest to the current joints)', ...
        nnz(Ipark.ok), Ipark.branch{:}), 0, iopNow, Vnow, 0, tproc, dt, a, camOff); %#ok<AGROW>
end
planOpt.seed = 11;
PLr = fp_arm_plan(QA(:,end), qPark, P, planOpt);
[QA, AI, T15] = planFrames(PLr, xPark, P);
st15 = PLr.stats;
for k = 1:size(QA,2)
    TL(end+1) = frame(P.q_home, QA(:,k), 15, sprintf(['Return path (RRT-Connect): speed %5.1f mm/s (SSM limit %5.1f)  ' ...
        'clearance %5.1f mm   remaining %5.1f mm'], AI(k).v*mm, PLr.vlim(k)*mm, PLr.clr(k), AI(k).d*mm), ...
        0, iopNow, Vnow, 0, tproc + (k-1)*dt, dt, AI(k), camOff); %#ok<AGROW>
end
tproc = tproc + T15; pathBack = AI(1).path; qEnd = QA(:,end); rawBack = rawTrace(PLr, P);
Tend = fp_arm_fk(qEnd, P); Tstart = fp_arm_fk(qPark, P);
retP = norm(Tend(1:3,4) - Tstart(1:3,4)); dRr = Tend(1:3,1:3).'*Tstart(1:3,1:3);
retR = acos(min(1, (trace(dRr) - 1)/2)); retQ = norm(qEnd - qPark);
parkedArm = armInfo(qEnd, xPark, P, 0, Ipark); parkedArm.br = 'parked';
for k = 1:round(2.6/dt)
    TL(end+1) = frame(P.q_home, qEnd, 15, sprintf(['Back at the start pose:  position error %.1e mm,  orientation error %.1e deg,  ' ...
        '|q_{end} - q_{start}| = %.1e rad'], retP*mm, rad2deg(retR), retQ), 0, iopNow, Vnow, 0, tproc, dt, parkedArm, camOff); %#ok<AGROW>
end
fprintf('Camera: entry-point error %.1f um | arm: docking residual %.1e mm, return error %.1e mm / %.1e deg, joint return %.1e rad\n', ...
    REG.err_entry*um, dockErrEst*mm, retP*mm, rad2deg(retR), retQ);
fprintf(['Arm planning (RRT-Connect, d_safe %g mm): out %d iter / %.1f s, min clearance %.1f mm, min stop margin %.1f mm, self %.1f mm | ' ...
    'back %d iter / %.1f s, min clearance %.1f mm, min stop margin %.1f mm, self %.1f mm\n'], planOpt.dSafe, st4.iter, st4.time, st4.minClr, ...
    st4.minSSM, st4.minSelf, st15.iter, st15.time, st15.minClr, st15.minSSM, st15.minSelf);
allQA = [TL.qa]; dQ = max(abs(diff(allQA, 1, 2)), [], 2)/dt;
fprintf('Arm peak joint speeds (deg/s): %s\n', sprintf('%.1f ', rad2deg(dQ)));
tdisp = cumsum([TL.pace]);                % display timeline (s)
if dry, assignin('base', 'TLdemo', TL); return; end

%% ---------------- 2. Figure layout (dark operating-room theme) ---------
bg = [0.05 0.07 0.10]; pan = [0.09 0.11 0.15]; fg = [0.88 0.91 0.95]; acc = [0.35 0.80 1.0]; gold = [1 0.85 0.3];
fig = figure('Name', 'MTR 5100 - camera + 4-DOF arm + 3-DOF RCM injection robot (live demo)', 'Color', bg, ...
    'Position', [30 40 1500 860], 'Visible', 'on', 'InvertHardcopy', 'off');
if ~usejava('desktop'), fig.Visible = 'off'; end
% every 3-D view lives in its own panel -> clipped to its own rectangle
pM = uipanel(fig, 'Position', [0.005 0.075 0.6 0.85], 'BorderType', 'none', 'BackgroundColor', bg);
pZ = uipanel(fig, 'Position', [0.61 0.47 0.22 0.45], 'BorderType', 'line', 'HighlightColor', [0.2 0.3 0.4], ...
    'BackgroundColor', bg);
axBG = axes(pM, 'Position', [0 0 1 1]);
gy = linspace(0, 1, 200).'; gx = linspace(-1, 1, 300);
G = 0.04 + 0.10*(1 - gy).^1.5 .* (1 - 0.6*gx.^2);
image(axBG, cat(3, 0.55*G, 0.75*G, 1.0*G)); axis(axBG, 'off');
% main cinematic view
axM = axes(pM, 'Position', [0 0 1 1], 'Color', 'none');
H = fp_scene(axM, P);
axis(axM, 'off'); daspect(axM, [1 1 1]);
xlim(axM, [-1300 1100]); ylim(axM, [-1300 1100]); zlim(axM, [-570 900]);
light(axM, 'Position', [0.4 -1 2], 'Style', 'infinite', 'Color', [1 0.97 0.92]);
light(axM, 'Position', [-1 0.8 0.5], 'Style', 'infinite', 'Color', [0.45 0.5 0.6]);
light(axM, 'Position', [0.2 1 -0.3], 'Style', 'infinite', 'Color', [0.25 0.35 0.5]);     % rim light
lighting(axM, 'gouraud'); axM.AmbientLightColor = [0.85 0.88 0.95];
hudAxis = plot3(axM, nan, nan, nan, '--', 'Color', [0.3 0.95 1], 'LineWidth', 1.3);
hudTgt  = plot3(axM, nan, nan, nan, 'o', 'Color', [0.3 0.95 1], 'MarkerSize', 7, 'LineWidth', 1.5);
% camera overlays: rays to the detected points, triangulated points, fitted limbus, entry estimate
Xm = REG.X*mm; cc = [REG.cams(1).Cc REG.cams(2).Cc]*mm;
rayIdx = 1:3:nP;
rayX = nan(1, 3*2*numel(rayIdx)); rayY = rayX; rayZ = rayX; m = 0;
for i = rayIdx
    for c = 1:2
        rayX(m+(1:3)) = [cc(1,c) Xm(1,i) nan]; rayY(m+(1:3)) = [cc(2,c) Xm(2,i) nan]; rayZ(m+(1:3)) = [cc(3,c) Xm(3,i) nan];
        m = m + 3;
    end
end
hRays = plot3(axM, nan, nan, nan, '-', 'Color', [0.3 1 0.5 0.35], 'LineWidth', 0.8);
hPts  = plot3(axM, nan, nan, nan, '.', 'Color', [0.3 1 0.5], 'MarkerSize', 9);
fc = REG.fitCircle3D*mm;
hFit  = plot3(axM, fc(1,:), fc(2,:), fc(3,:), '-', 'Color', [0.3 0.95 1], 'LineWidth', 2, 'Visible', 'off');
pEst = REG.entry*mm;
hEntry = plot3(axM, pEst(1), pEst(2), pEst(3), 's', 'Color', gold, 'MarkerFaceColor', gold, 'MarkerSize', 6, 'Visible', 'off');
hEntryT = text(axM, pEst(1) + 3, pEst(2), pEst(3) + 3, 'entry (camera)', 'Color', gold, 'FontSize', 9, ...
    'FontWeight', 'bold', 'Visible', 'off');
% arm overlays: planned Cartesian path, executed trail, frame triads, start marker
pathPlan  = plot3(axM, nan, nan, nan, ':', 'Color', gold, 'LineWidth', 1.8);
pathRaw   = plot3(axM, nan, nan, nan, '-', 'Color', [1 0.35 0.3 0.55], 'LineWidth', 1.2);   % RRT path before shortcut + smoothing
pathTrail = animatedline(axM, 'Color', [0.3 0.95 1], 'LineWidth', 2.4);
triCur  = triad(axM, 70, 2.2, '');
triGoal = triad(axM, 70, 1.4, 'IK target');
Tstart_mm = Tstart; pS = Tstart_mm(1:3,4).'*mm;
plot3(axM, pS(1), pS(2), pS(3), 'd', 'MarkerSize', 9, 'LineWidth', 1.6, 'Color', gold, 'MarkerFaceColor', gold*0.4);
text(axM, pS(1), pS(2), pS(3) - 45, 'START', 'Color', gold, 'FontWeight', 'bold', 'FontSize', 10, ...
    'HorizontalAlignment', 'center');
% camera keyframes per step: {target, az, el, distance, view angle}
pEm = P.eye.pE.'*mm;
wide = {[80 -100 220], 60, 20, 2800, 26};
KF = {{[-60 60 170], 35, 16, 2800, 22}, ...        % 1 camera + eye
      {[0 0 8], 25, 55, 2800, 1.5}, ...            % 2 eye close-up with the fitted limbus
      wide, ...                                    % 3 IK + ghost
      {[60 -80 220], 85, 22, 2800, 24}, ...        % 4 transport
      {[-25 -10 60], 155, 28, 2800, 9.5}, ...      % 5 docking (head side, column out of the way)
      {[-10 -5 30], 140, 30, 1400, 11}, ...        % 6 align
      {[3 0 12],    128, 30, 1400, 6.5}, ...       % 7 approach
      {[3 0 12],    124, 30, 1400, 6.5}, ...       % 8 insertion close-up
      {[3 0 12],    132, 33, 1400, 7.0}, ...       % 9 injection A
      {[3 0 12],    150, 40, 1400, 7.5}, ...       % 10 re-orientation (slow orbit)
      {[3 0 12],    162, 40, 1400, 7.0}, ...       % 11 injection B (surgeon's view)
      {pEm + [-2 0 2], 188, 40, 1400, 5.5}, ...    % 12 twitch (tight)
      {[-25 -25 40], 215, 32, 1400, 12}, ...       % 13 retract
      {[-20 -20 70], 245, 26, 2800, 7.5}, ...      % 14 undock
      {[80 -100 220], 300, 22, 2800, 26}};         % 15 IK return (wide)
fp_set_camera(axM, KF{1}{:});
% right-top panel: camera image (steps 1-2) / eye x-ray close-up (later)
pZc = uipanel(pZ, 'Position', [0 0 1 1], 'BorderType', 'none', 'BackgroundColor', bg);
pZx = uipanel(pZ, 'Position', [0 0 1 1], 'BorderType', 'none', 'BackgroundColor', bg, 'Visible', 'off');
axC = axes(pZc, 'Position', [0.04 0.04 0.92 0.84], 'Color', [0.12 0.12 0.13], 'XColor', 'none', 'YColor', 'none');
hold(axC, 'on'); axis(axC, 'ij'); axis(axC, 'equal');
ctr = REG.img.centre; win = 165;
xlim(axC, ctr(1) + [-win win]); ylim(axC, ctr(2) + [-win win]*0.95);
[camImg, pp, hw] = fp_camera_image(P, REG, qPark);      % what the left camera actually sees (rendered 3-D scene)
image(axC, 'XData', pp(1) + [-hw hw], 'YData', pp(2) + [-hw hw], 'CData', camImg);
cDet = plot(axC, nan, nan, '+', 'Color', [0.2 1 0.4], 'MarkerSize', 7, 'LineWidth', 1.3);
cFit = plot(axC, REG.img.fit(1,:), REG.img.fit(2,:), '-', 'Color', [0.3 0.95 1], 'LineWidth', 1.8, 'Visible', 'off');
cEnt = plot(axC, REG.img.entry(1), REG.img.entry(2), 's', 'Color', gold, 'MarkerFaceColor', gold, 'MarkerSize', 7, 'Visible', 'off');
cTxt = text(axC, ctr(1) - win + 8, ctr(2) - 0.95*win + 8, '', 'Color', fg, 'FontName', 'Consolas', 'FontSize', 8.5, ...
    'VerticalAlignment', 'top', 'BackgroundColor', [0.05 0.07 0.10], 'Margin', 3);
axTc = axes(pZc, 'Position', [0 0.9 1 0.1], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
text(axTc, 0.5, 0.5, sprintf('LEFT CAMERA IMAGE  (%d x %d px)', P.cam.res), 'Color', fg, 'FontSize', 10, ...
    'FontWeight', 'bold', 'HorizontalAlignment', 'center');
axZ = axes(pZx, 'Position', [0 0 1 0.92], 'Color', bg);
HZ = fp_scene(axZ, P, struct('head', false, 'setup', false, 'xray', true));
axis(axZ, 'equal'); axis(axZ, 'off'); view(axZ, -28, 14);
xlim(axZ, [-13 13]); ylim(axZ, [-13 13]); zlim(axZ, [-13 15]);
light(axZ, 'Position', [-1 -0.5 1], 'Style', 'infinite');
light(axZ, 'Position', [1 1 0.3], 'Style', 'infinite', 'Color', [0.4 0.45 0.55]);
lighting(axZ, 'gouraud');
axT = axes(pZx, 'Position', [0 0.92 1 0.08], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1]);
text(axT, 0.5, 0.5, 'EYE CLOSE-UP  (x-ray view)', 'Color', fg, 'FontSize', 10, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'center');
tipTrail = animatedline(axZ, 'Color', acc, 'LineWidth', 1.8);
% joint readout panel
axJ = axes(fig, 'Position', [0.84 0.47 0.15 0.45], 'Color', pan, 'XColor', [0.2 0.3 0.4], 'YColor', [0.2 0.3 0.4], ...
    'XLim', [0 1], 'YLim', [0 1], 'Box', 'on', 'XTick', [], 'YTick', []);
txtJ = text(axJ, 0.06, 0.97, '', 'FontName', 'Consolas', 'FontSize', 9.2, 'VerticalAlignment', 'top', 'Color', fg);
% live plots
axF = axes(fig, 'Position', [0.645 0.265 0.325 0.15]); hold(axF, 'on');
yyaxis(axF, 'left'); lnD = animatedline(axF, 'Color', acc, 'LineWidth', 1.6); ylabel(axF, 'depth q_3 (mm)');
ylim(axF, [-16 10]);
yyaxis(axF, 'right'); lnF = animatedline(axF, 'Color', [1 0.55 0.2], 'LineWidth', 1.3); ylabel(axF, 'force (N)');
ylim(axF, [0 0.25]); xlim(axF, [0 tdisp(end)]); title(axF, 'Insertion depth and needle force', 'Color', fg);
axP = axes(fig, 'Position', [0.645 0.06 0.325 0.15]); hold(axP, 'on');
yyaxis(axP, 'left'); lnP = animatedline(axP, 'Color', [0.95 0.45 0.9], 'LineWidth', 1.6); ylabel(axP, 'IOP (mmHg)');
ylim(axP, [10 30]); yline(axP, P.iop.Plim, '--', 'Color', [1 0.3 0.3]);
yyaxis(axP, 'right'); lnV = animatedline(axP, 'Color', [0.4 0.9 0.5], 'LineWidth', 1.6); ylabel(axP, 'volume (\muL)');
ylim(axP, [0 55]); xlim(axP, [0 tdisp(end)]); xlabel(axP, 'demo timeline (s)  -  injections shown in time-lapse');
title(axP, 'Intraocular pressure and delivered dose', 'Color', fg);
for a = [axF axP]
    set(a, 'Color', pan, 'XColor', fg, 'GridColor', [0.6 0.7 0.8], 'GridAlpha', 0.15, 'FontSize', 9);
    grid(a, 'on'); yyaxis(a, 'left'); a.YColor = fg; yyaxis(a, 'right'); a.YColor = fg;
end
% stage strip (camera + arm / RCM robot / arm)
t1 = tdisp(find([TL.ph] == 6, 1)) - dt; t2 = tdisp(find([TL.ph] == 14, 1)) - dt; tc = tdisp(find([TL.ph] == 3, 1)) - dt;
for a = [axF axP]
    yyaxis(a, 'left'); yl = ylim(a);
    patch(a, [0 t1 t1 0], yl([1 1 2 2]), gold, 'FaceAlpha', 0.07, 'EdgeColor', 'none');
    patch(a, [t2 tdisp(end) tdisp(end) t2], yl([1 1 2 2]), gold, 'FaceAlpha', 0.07, 'EdgeColor', 'none');
end
text(axF, 0.5*tc, -13, 'CAMERA', 'Color', [0.3 1 0.5], 'FontSize', 8, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(axF, 0.5*(tc + t1), -13, 'ARM', 'Color', gold, 'FontSize', 8, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(axF, 0.5*(t2 + tdisp(end)), -13, 'ARM', 'Color', gold, 'FontSize', 8, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(axF, 0.5*(t1 + t2), -13, 'RCM ROBOT', 'Color', acc, 'FontSize', 8, 'HorizontalAlignment', 'center', 'FontWeight', 'bold');
% banners
ttl = annotation(fig, 'textbox', [0.01 0.93 0.6 0.06], 'String', '', 'FontSize', 16, 'FontWeight', 'bold', ...
    'EdgeColor', 'none', 'Color', acc);
sub = annotation(fig, 'textbox', [0.01 0.012 0.6 0.055], 'String', '', 'FontSize', 12, 'EdgeColor', 'none', ...
    'Color', [1 0.8 0.35]);
annotation(fig, 'textbox', [0.62 0.93 0.37 0.06], 'String', 'MTR 5100  |  Camera + 4-DOF arm + 3-DOF RCM Injection Robot', ...
    'FontSize', 10.5, 'EdgeColor', 'none', 'HorizontalAlignment', 'right', 'Color', [0.6 0.66 0.75]);

if record
    vw = VideoWriter(fullfile(here, '..', 'results', 'demo_procedure.mp4'), 'MPEG-4');
    vw.FrameRate = 25; vw.Quality = 92; open(vw);
end

%% ---------------- 3. Play -----------------------------------------------
fprintf('Playing %d frames...\n', numel(TL));
phs = [TL.ph]; snapIdx = zeros(1,nPh); phStart = zeros(1,nPh); phLen = zeros(1,nPh);
for p = 1:nPh
    ii = find(phs == p); snapIdx(p) = ii(max(1, round(0.8*numel(ii))));
    phStart(p) = ii(1); phLen(p) = numel(ii);
end
ii = find(phs == 12); snapIdx(12) = ii(round(0.45*numel(ii)));
ii = find(phs == 4);  snapIdx(4)  = ii(round(0.5*numel(ii)));
ii = find(phs == 15); snapIdx(15) = ii(round(0.55*numel(ii)));
lastPh = 0; camFrom = KF{1};
tgtA = (P.T_E0*[fp_fk_geom(qA); 1]).'*mm; tgtB = (P.T_E0*[fp_fk_geom(qB); 1]).'*mm;
[~, TgD] = fp_arm_fk(qDockSol, P); [~, TgP] = fp_arm_fk(qPark, P);
TFgoalD = fp_arm_fk(qDockSol, P);
for k = 1:numel(TL)
    t0 = tic;
    f = TL(k);
    if ~isvalid(fig), return; end
    TE0now = fp_scene_update(H,  f.q, P, f.th, f.V, f.qa);
    fp_scene_update(HZ, f.q, P, f.th, f.V, f.qa);
    if f.ph ~= lastPh
        ttl.String = ['Step ' phaseName{f.ph}];
        if lastPh > 0, camFrom = camNow; end
        lastPh = f.ph;
        if f.ph == 12 || f.ph == 13, clearpoints(tipTrail); end
        pZc.Visible = matlab.lang.OnOffSwitchState(f.ph <= 2);
        pZx.Visible = matlab.lang.OnOffSwitchState(f.ph > 2);
        clearpoints(pathTrail); setPath(pathRaw, []);
        switch f.ph
            case {3, 4}, setPath(pathPlan, pathOut); showGhost(H.ghost, TgD, true); setTriad(triGoal, TFgoalD);
                         setPath(pathRaw, rawOut);
            case 5,      setPath(pathPlan, pathDock); showGhost(H.ghost, TgD, false); setTriad(triGoal, TFgoalD);
            case 14,     setPath(pathPlan, flip(pathDock, 2)); showGhost(H.ghost, TgD, false); setTriad(triGoal, []);
            case 15,     setPath(pathPlan, pathBack); showGhost(H.ghost, TgP, true); setTriad(triGoal, Tstart);
                         setPath(pathRaw, rawBack);
            otherwise,   setPath(pathPlan, []); showGhost(H.ghost, TgD, false); setTriad(triGoal, []);
        end
    end
    % ---- camera overlays ----
    nS = f.cam.n;
    if f.cam.rays && nS > 0
        nr = nnz(rayIdx <= nS)*6;
        set(hRays, 'XData', rayX(1:nr), 'YData', rayY(1:nr), 'ZData', rayZ(1:nr));
    else
        set(hRays, 'XData', nan, 'YData', nan, 'ZData', nan);
    end
    if nS > 0
        set(hPts, 'XData', Xm(1,1:nS), 'YData', Xm(2,1:nS), 'ZData', Xm(3,1:nS));
        set(cDet, 'XData', REG.uv{1}(1,1:nS), 'YData', REG.uv{1}(2,1:nS));
    elseif ~f.cam.fit
        set(hPts, 'XData', nan, 'YData', nan, 'ZData', nan);
    end
    set([hFit cFit cEnt], 'Visible', matlab.lang.OnOffSwitchState(f.cam.fit && f.ph <= 2));
    set([hEntry hEntryT], 'Visible', matlab.lang.OnOffSwitchState(f.cam.fit));
    if f.ph <= 2
        if f.cam.fit
            cTxt.String = sprintf('limbus fit r = %.3f mm\nentry error  = %.0f um\naxis error   = %.2f deg', ...
                REG.limbus_r*mm, REG.err_entry*um, rad2deg(REG.err_axis));
        else
            cTxt.String = sprintf('edge points: %d / %d\nnoise %.2f px', nS, nP, P.cam.sd_px);
        end
    end
    % ---- cinematic camera: ease from the previous shot to this step's shot ----
    nT = max(1, min(45, round(0.6*phLen(f.ph))));
    s = min(1, (k - phStart(f.ph))/nT); s = s*s*(3 - 2*s);
    kf = KF{f.ph};
    camNow = {camFrom{1} + s*(kf{1} - camFrom{1}), camFrom{2} + s*(kf{2} - camFrom{2}), ...
              camFrom{3} + s*(kf{3} - camFrom{3}), camFrom{4} + s*(kf{4} - camFrom{4}), camFrom{5}*(kf{5}/camFrom{5})^s};
    fp_set_camera(axM, camNow{:});
    % ---- arm overlays: flange triad and executed path ----
    TF = fp_arm_fk(f.qa, P);
    setTriad(triCur, TF);
    if any(f.ph == [4 5 14 15]), addpoints(pathTrail, TF(1,4)*mm, TF(2,4)*mm, TF(3,4)*mm); end
    % ---- HUD: planned needle axis and target ----
    [~, u0] = fp_fk_geom(f.q); uE = TE0now(1:3,1:3)*u0; pRCM = TE0now(1:3,4).'*mm;
    if f.ph >= 6 && f.ph <= 11
        tg = tgtA; if f.ph >= 10, tg = tgtB; end
        a = pRCM - 40*uE.'; b = pRCM + 7*uE.';
        set(hudAxis, 'XData', [a(1) b(1)], 'YData', [a(2) b(2)], 'ZData', [a(3) b(3)]);
        set(hudTgt, 'XData', tg(1), 'YData', tg(2), 'ZData', tg(3));
    else
        set(hudAxis, 'XData', nan, 'YData', nan, 'ZData', nan); set(hudTgt, 'XData', nan, 'YData', nan, 'ZData', nan);
    end
    sub.String = f.msg;
    pt = TE0now*[fp_fk_geom(f.q); 1];
    if f.q(3) > 0 && f.ph >= 6 && f.ph <= 13, addpoints(tipTrail, pt(1)*mm, pt(2)*mm, pt(3)*mm); end
    addpoints(lnD, tdisp(k), f.q(3)*mm); addpoints(lnF, tdisp(k), f.F);
    addpoints(lnP, tdisp(k), f.iop);    addpoints(lnV, tdisp(k), f.V);
    qa = rad2deg(f.qa);
    txtJ.String = sprintf(['CAMERA REGISTRATION\n entry error %5.0f um\n\n4-DOF ARM  (%s)\n q1 base     %7.1f deg\n' ...
        ' q2 shoulder %7.1f deg\n q3 elbow    %7.1f deg\n q4 wrist    %7.1f deg\n IK residual %.0e mm\n |det J|     %.2e\n\n' ...
        '3-DOF RCM ROBOT\n q1 yaw      %7.2f deg\n q2 pitch    %7.2f deg\n q3 depth    %7.2f mm\n\nTIME  t = %6.1f s'], ...
        REG.err_entry*um, f.arm.br, qa, max(f.arm.res*mm, 1e-16), f.arm.w, rad2deg(f.q(1)), rad2deg(f.q(2)), f.q(3)*mm, f.tp);
    drawnow limitrate;
    if record, writeVideo(vw, getframe(fig)); end
    if snaps && any(k == snapIdx)
        drawnow;
        ws = warning('off', 'all');
        print(fig, fullfile(here, '..', 'figures', sprintf('D%02d_demo_step%d', f.ph, f.ph)), '-dpng', '-r110');
        warning(ws);
    end
    if ~fast && ~record, pause(max(0, f.pace - toc(t0))); end
end
drawnow;
if record, close(vw); fprintf('Video saved to results/demo_procedure.mp4\n'); end
fprintf('Demo finished.\n');
end

%% ======================= timeline helpers ==============================
function f = frame(q, qa, ph, msg, F, iop, V, th, tp, pace, arm, cam)
f = struct('q', q, 'qa', qa, 'ph', ph, 'msg', msg, 'F', F, 'iop', iop, 'V', V, 'th', th, 'tp', tp, ...
    'pace', pace, 'arm', arm, 'cam', cam);
end

function x = taskOf(T, P)
% task coordinates of a flange pose: [position in {E}; pitch angle phi]
Tb = P.arm.T_Eb\T; pb = Tb(1:3,4); xb = Tb(1:3,1);
h = [pb(1); pb(2); 0]; h = h/norm(h);
x = [T(1:3,4); atan2(xb(3), h.'*xb)];
end

function x = taskForRCM(pR, phi, P)
% flange task pose [p4; phi] whose RCM point (dF along the flange x-axis)
% lands exactly on pR: the flange x-axis direction depends on the arm-plane
% direction (base yaw), so solve the small fixed-point problem
%   p4 = pR - dF * (cos(phi) h(p4) + sin(phi) up)
Tb = P.arm.T_Eb; up = Tb(1:3,3); p4 = pR;
for it = 1:20
    pb = Tb\[p4; 1]; h = Tb(1:3,1:3)*[pb(1); pb(2); 0]; h = h/norm(h);
    p4 = pR - P.arm.dF*(cos(phi)*h + sin(phi)*up);
end
x = [p4; phi];
end

function [QA, AI, T] = armPath(x0, x1, vias, qseed, P, vmax, dt)
% Cartesian path in the arm's task space [p; phi]: the flange position
% follows a pchip spline through the via points (chord-length parameter),
% phi changes linearly; both follow one quintic time law whose peak speed
% (1.875 L/T) equals vmax. Closed-form IK at every frame, closest branch.
pts = [x0(1:3) vias x1(1:3)];
c = [0 cumsum(vecnorm(diff(pts, 1, 2)))];
sd = linspace(0, c(end), 400);
dense = [pchip(c, pts(1,:), sd); pchip(c, pts(2,:), sd); pchip(c, pts(3,:), sd)];
cs = [0 cumsum(vecnorm(diff(dense, 1, 2)))]; L = cs(end);
T = max(1.875*L/vmax, 1.5); t = 0:dt:T; if t(end) < T, t = [t T]; end
s = 10*(t/T).^3 - 15*(t/T).^4 + 6*(t/T).^5;
QA = zeros(4, numel(t)); q = qseed; pPrev = x0(1:3);
AI = repmat(struct('res', 0, 'w', 0, 'v', 0, 'br', '', 'd', 0, 'path', []), 1, numel(t));
for k = 1:numel(t)
    p = [interp1(cs, dense(1,:), s(k)*L); interp1(cs, dense(2,:), s(k)*L); interp1(cs, dense(3,:), s(k)*L)];
    xk = [p; x0(4) + s(k)*(x1(4) - x0(4))];
    [q, I] = fp_arm_ik(xk, P, q);
    QA(:,k) = q;
    v = norm(p - pPrev)/dt; pPrev = p;
    AI(k) = armInfo(q, xk, P, v, I);
    AI(k).d = norm(x1(1:3) - p);
end
AI(1).path = dense;
end

function [QA, AI, T] = planFrames(PL, xGoal, P)
% timeline frames of a planned arm trajectory (FP_ARM_PLAN)
QA = PL.Q; T = PL.T; n = size(QA,2);
AI = repmat(struct('res', 0, 'w', 0, 'v', 0, 'br', '', 'd', 0, 'path', []), 1, n);
for k = 1:n
    TF = fp_arm_fk(QA(:,k), P);
    a = armInfo(QA(:,k), [TF(1:3,4); 0], P, max(PL.v(k), PL.vR(k)), struct('branch', {{'planned', 'RRT-Connect'}}));
    a.d = norm(xGoal(1:3) - TF(1:3,4)); AI(k) = a;
end
AI(1).path = PL.pF;
end

function pts = rawTrace(PL, P)
% flange trace of the raw RRT-Connect path (straight joint-space edges)
R = PL.raw; pts = zeros(3,0);
for i = 1:size(R,2)-1
    for s = linspace(0, 1, 12)
        TF = fp_arm_fk(R(:,i) + s*(R(:,i+1) - R(:,i)), P); pts(:,end+1) = TF(1:3,4); %#ok<AGROW>
    end
end
end

function a = armInfo(q, xg, P, v, I)
% IK residual, task-Jacobian manipulability |det J| and branch label for the HUD
[TF, Ta] = fp_arm_fk(q, P);
o = [P.arm.T_Eb(1:3,4) squeeze(Ta(1:3,4,1:3))]; z = [P.arm.T_Eb(1:3,3) squeeze(Ta(1:3,3,1:3))];
J = zeros(4);
for i = 1:4, J(1:3,i) = cross(z(:,i), TF(1:3,4) - o(:,i)); end
J(4,:) = [0 1 1 1];                                   % phi = q2 + q3 + q4
res = norm(TF(1:3,4) - xg(1:3));
a = struct('res', res, 'w', abs(det(J)), 'v', v, 'br', strjoin(I.branch, ', '), 'd', 0, 'path', []);
end

%% ======================= overlay helpers ===============================
function h = triad(ax, len, lw, label)
h.T = hgtransform('Parent', ax, 'Visible', 'off');
cols = [1 0.3 0.3; 0.3 1 0.4; 0.35 0.6 1];
for i = 1:3
    e = zeros(3,1); e(i) = len;
    plot3([0 e(1)], [0 e(2)], [0 e(3)], '-', 'Color', cols(i,:), 'LineWidth', lw, 'Parent', h.T);
end
if ~isempty(label)
    text(0, 0, -25, label, 'Color', [0.35 0.8 1], 'FontSize', 9, 'FontWeight', 'bold', 'Parent', h.T);
end
end

function setTriad(h, T)
if isempty(T), h.T.Visible = 'off'; return; end
h.T.Matrix = [T(1:3,1:3) T(1:3,4)*1e3; 0 0 0 1]; h.T.Visible = 'on';
end

function setPath(h, pts)
if isempty(pts), set(h, 'XData', nan, 'YData', nan, 'ZData', nan); return; end
set(h, 'XData', pts(1,:)*1e3, 'YData', pts(2,:)*1e3, 'ZData', pts(3,:)*1e3);
end

function showGhost(G, Ta, on)
if isempty(G), return; end
for i = 1:4, G.T(i+1).Matrix = [Ta(1:3,1:3,i) Ta(1:3,4,i)*1e3; 0 0 0 1]; end
set(G.T(2:end), 'Visible', matlab.lang.OnOffSwitchState(on));
end
