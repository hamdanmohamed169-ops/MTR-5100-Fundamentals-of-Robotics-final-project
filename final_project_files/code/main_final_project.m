%% MTR 5100 - Fundamentals of Robotics - FINAL PROJECT
%  Design, Kinematic Modelling and Trajectory Control of a Micro-Precision
%  Surgical Robot with Remote Centre of Motion (RCM) for Intravitreal
%  Retinal Injections:
%    - stereo camera that locates the eye (registration)
%    - 4-DOF positioning arm that carries the injection unit to the eye
%    - 3-DOF RCM injection robot (yaw, pitch, insertion)
%
%  Author : Mohamed Abdelaziz Hamdan   (ID 481019296)
%  MATLAB : R2024a, base MATLAB only (no toolboxes)
%
%  HOW TO RUN: open this folder in MATLAB and run  main_final_project
%  All figures are written to ../figures and the console log to
%  ../results/run_log.txt. Runtime is about 11-26 minutes (machine dependent).
%
%  Course coverage (Module -> section of this script)
%   M1 Robot anatomy / workspace ............. S1
%   M2 Actuation, transmissions, end effector  S9, S9b (arm drives, brakes) + reflected inertia
%   M3 Sensing (camera, encoders, force, pressure, eye tracker) S2, S5, S8, S11, S12
%   M4 Frames and homogeneous transforms, DH . S2, S3
%   M5 Forward / inverse kinematics (RCM robot + 4-DOF arm) S3
%   M6 Jacobian, singularities, statics, dynamics S4, S9b (arm dynamics)
%   M7 Trajectories and control .............. S6, S8
%   M8 Motion planning (C-space, A*, RRT, smoothing) S7, S7b (arm: RRT-Connect around the patient)
%   M9 Safety, ethics, professional practice . S10, S12, S13 (risk register,
%      stopping distance, ethics), S14 (summary, FMEA)

clear; clc; close all;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
cd(here);
figDir = fullfile(here, '..', 'figures');  resDir = fullfile(here, '..', 'results');
if ~exist(figDir,'dir'), mkdir(figDir); end
if ~exist(resDir,'dir'), mkdir(resDir); end
logFile = fullfile(resDir, 'run_log.txt');
if exist(logFile,'file'), delete(logFile); end
diary(logFile);
tStart = tic;

% ---- Live 3-D prototype demonstration of the full procedure ----------
% Set RUN_DEMO = false to skip the animation and run only the analysis.
RUN_DEMO = true;
if RUN_DEMO && usejava('desktop')
    demo_live();
    fprintf('\nDemo done - starting the full engineering analysis...\n\n');
end

P = fp_params();
rng(P.seed);
mm = 1e3; um = 1e6;                       % display conversions
R0E = P.T_E0(1:3,1:3); oE = P.T_E0(1:3,4);
toE = @(p0) oE + R0E*p0;                  % {0} -> {E}
to0 = @(pE) R0E.'*(pE - oE);              % {E} -> {0}
SUM = {};                                 % summary rows {metric, target, achieved, pass}

fprintf('=====================================================================\n');
fprintf(' MTR 5100 Final Project - camera + 4-DOF arm + 3-DOF RCM injection robot\n');
fprintf('=====================================================================\n\n');

%% ======================================================================
%  S1. SYSTEM DESCRIPTION AND SAFE WORKSPACE  (Module 1)
%  ======================================================================
fprintf('--- S1. Robot anatomy and workspace ---\n');
fprintf('Architecture: R(yaw)-R(pitch)-P(insertion), 3 DOF spherical robot, RCM at pars plana\n');
fprintf('Entry point (RCM) in {E}: [%.2f %.2f %.2f] mm, %.1f deg from the anterior pole\n', ...
    oE*mm, rad2deg(P.eye.theta_entry));
fprintf('Joint limits: q1 [%g %g] deg, q2 [%g %g] deg, q3 [%g %g] mm\n', ...
    rad2deg(P.qmin(1)), rad2deg(P.qmax(1)), rad2deg(P.qmin(2)), rad2deg(P.qmax(2)), ...
    P.qmin(3)*mm, P.qmax(3)*mm);

% Cartesian grid over the vitreous: which points can the tip reach SAFELY?
h = 0.5e-3;  g = -P.eye.Rret:h:P.eye.Rret;
[X,Y,Z] = ndgrid(g,g,g);
inside = X.^2 + Y.^2 + Z.^2 <= P.eye.Rret^2;
lensIn = ((X-P.eye.lens_c(1))/P.eye.lens_ax(1)).^2 + ((Y-P.eye.lens_c(2))/P.eye.lens_ax(2)).^2 + ...
         ((Z-P.eye.lens_c(3))/P.eye.lens_ax(3)).^2 <= 1;
vit = inside & ~lensIn;                       % vitreous cavity (approx.)
idxV = find(vit);
reach = false(numel(idxV),1); safe = reach;
Pw = P; Pw.eye.keepout_r = -1;               % anatomical workspace: patient keep-out zone not applied here
for k = 1:numel(idxV)
    pe = [X(idxV(k)); Y(idxV(k)); Z(idxV(k))];
    [qk, ik] = fp_ik(to0(pe), P, P.q_home);
    if isempty(qk), continue; end
    reach(k) = ik.err(ik.selected) < 1e-12;
    safe(k)  = reach(k) && qk(3) <= P.depth_soft && fp_lens_check(qk, Pw);
end
Vvit = numel(idxV)*h^3; Vreach = sum(reach)*h^3; Vsafe = sum(safe)*h^3;
fprintf('Vitreous volume (model)            : %7.0f mm^3\n', Vvit*1e9);
fprintf('Kinematically reachable volume     : %7.0f mm^3 (%.1f %%)\n', Vreach*1e9, 100*Vreach/Vvit);
fprintf('Safe reachable volume (clearances, depth limit): %7.0f mm^3 (%.1f %%)\n\n', Vsafe*1e9, 100*Vsafe/Vvit);

fig = figure('Color','w','Position',[50 50 1100 480]);
subplot(1,2,1); hold on; axis equal; grid on; box on;
pts = [X(idxV(safe)) Y(idxV(safe)) Z(idxV(safe))]*mm;
pts = pts(1:4:end,:);
scatter3(pts(:,1),pts(:,2),pts(:,3),3,pts(:,3),'filled','MarkerFaceAlpha',0.25);
[sx,sy,sz] = sphere(30);
surf(P.eye.R*sx*mm,P.eye.R*sy*mm,P.eye.R*sz*mm,'FaceAlpha',0.05,'EdgeAlpha',0.1,'FaceColor',[0.5 0.5 0.5]);
surf(P.eye.lens_ax(1)*sx*mm+P.eye.lens_c(1)*mm, P.eye.lens_ax(2)*sy*mm, ...
     P.eye.lens_ax(3)*sz*mm+P.eye.lens_c(3)*mm,'FaceAlpha',0.3,'EdgeColor','none','FaceColor',[0.9 0.6 0.2]);
plot3(oE(1)*mm,oE(2)*mm,oE(3)*mm,'rp','MarkerSize',14,'MarkerFaceColor','r');
xlabel('x_E (mm)'); ylabel('y_E (mm)'); zlabel('z_E (mm)'); view(35,20);
title('Safe reachable tip workspace inside the eye');
legend({'safe tip','globe','lens','RCM entry'},'Location','southoutside');
subplot(1,2,2); hold on; axis equal; grid on; box on;
sl = abs(Y(idxV)) < h/2;                          % slice yE = 0
plot(X(idxV(sl & ~reach))*mm, Z(idxV(sl & ~reach))*mm,'.','Color',[0.8 0.8 0.8]);
plot(X(idxV(sl & reach & ~safe))*mm, Z(idxV(sl & reach & ~safe))*mm,'.','Color',[0.95 0.6 0.5]);
plot(X(idxV(sl & safe))*mm, Z(idxV(sl & safe))*mm,'.','Color',[0.2 0.6 0.3]);
tt = linspace(0,2*pi,200);
plot(P.eye.R*cos(tt)*mm, P.eye.R*sin(tt)*mm,'k-');
plot(P.eye.Rret*cos(tt)*mm, P.eye.Rret*sin(tt)*mm,'k:');
plot(P.eye.lens_ax(1)*cos(tt)*mm, (P.eye.lens_c(3)+P.eye.lens_ax(3)*sin(tt))*mm,'-','Color',[0.9 0.5 0.1],'LineWidth',1.5);
plot(oE(1)*mm,oE(3)*mm,'rp','MarkerSize',12,'MarkerFaceColor','r');
xlabel('x_E (mm)'); ylabel('z_E (mm)');
title('Slice y_E = 0: unreachable / reachable-unsafe / safe');
legend({'not reachable','violates clearance','safe'},'Location','southoutside');
saveFig(fig, figDir, 'F01_workspace');

%% ======================================================================
%  S2. COORDINATE FRAMES AND TRANSFORMS  (Module 4)
%  ======================================================================
fprintf('--- S2. Frames {E} (eye) and {0} (robot base at RCM) ---\n');
Rchk = R0E.'*R0E;
Tinv = [R0E.' -R0E.'*oE; 0 0 0 1];
fprintf('||R''R - I|| = %.2e,  det(R) = %.12f,  ||T*T^-1 - I|| = %.2e\n', ...
    norm(Rchk-eye(3)), det(R0E), norm(P.T_E0*Tinv - eye(4)));
% clinical targets defined in the eye frame
qA = [deg2rad(-20); deg2rad(100); 7e-3];    % configuration used for site A
qB = [deg2rad( 20); deg2rad(100); 7e-3];    % configuration used for site B
tgtE = [toE(fp_fk_geom(qA)) toE(fp_fk_geom(qB)) [0;0;-1e-3]];
tgtName = {'Site A','Site B','Deep central (unsafe test)'};
fprintf('%-28s %-26s %-26s\n','Target','p in {E} (mm)','p in {0} (mm)');
for k = 1:3
    fprintf('%-28s [%6.2f %6.2f %6.2f]     [%6.2f %6.2f %6.2f]\n', tgtName{k}, tgtE(:,k)*mm, to0(tgtE(:,k))*mm);
end
fprintf('Round trip {E}->{0}->{E} error: %.2e m\n\n', norm(toE(to0(tgtE(:,1))) - tgtE(:,1)));
SUM(end+1,:) = {'Frame transform validity ||R''R-I||', '< 1e-12', sprintf('%.1e', norm(Rchk-eye(3))), norm(Rchk-eye(3))<1e-12};

% ---- Eye registration with the stereo camera (Modules 3-4) ----
fprintf('--- S2b. Locating the eye with the stereo camera ---\n');
fprintf('Cameras: f = %d px, %dx%d px, baseline %.0f mm, %.0f mm from the eye, pixel noise %.2f px, %d limbus points\n', ...
    P.cam.f, P.cam.res, P.cam.base*mm, norm(P.cam.C - P.cam.target)*mm, P.cam.sd_px, P.cam.nPts);
REG = fp_eye_register(P, P.seed);
fprintf('Nominal run: limbus radius error %.1f um, eye-centre error %.1f um, optical-axis error %.3f deg\n', ...
    REG.err_radius*um, REG.err_centre*um, rad2deg(REG.err_axis));
fprintf('             entry-point (RCM target) error %.1f um, robot-base frame orientation error %.3f deg\n', ...
    REG.err_entry*um, rad2deg(REG.err_frame));
nReg = 300; eReg = zeros(1,nReg); aReg = eReg;
for k = 1:nReg
    Rk = fp_eye_register(P, P.seed + 1000 + k); eReg(k) = Rk.err_entry; aReg(k) = Rk.err_frame;
end
fprintf('Monte Carlo (%d trials): entry-point error median %.1f um, 95th pct %.1f um, max %.1f um (requirement %.0f um)\n\n', ...
    nReg, median(eReg)*um, prctile95(eReg)*um, max(eReg)*um, P.cam.req*um);
SUM(end+1,:) = {'Eye registration: entry-point error (95th pct)', sprintf('<= %.2f mm', P.cam.req*mm), ...
    sprintf('%.3f mm', prctile95(eReg)*mm), prctile95(eReg) <= P.cam.req};
fig = figure('Color','w','Position',[50 50 1150 450]);
subplot(1,2,1); hold on; axis ij; axis equal; box on;
qParkCam = fp_arm_ik([0.33; -0.25; 0.40; -pi/2], P, deg2rad([60; 80; -100; -70]));   % arm parked (out of view)
[camImg, ppC, hwC] = fp_camera_image(P, REG, qParkCam);       % rendered frame of the left camera
image(ppC(1) + [-hwC hwC], ppC(2) + [-hwC hwC], camImg);
plot(REG.uv{1}(1,:), REG.uv{1}(2,:), 'g+', 'MarkerSize', 7, 'LineWidth', 1.2);
plot(REG.img.fit(1,:), REG.img.fit(2,:), 'c-', 'LineWidth', 1.5);
plot(REG.img.entry(1), REG.img.entry(2), 'ys', 'MarkerFaceColor', 'y', 'MarkerSize', 8);
plot(REG.img.entryTrue(1), REG.img.entryTrue(2), 'r+', 'MarkerSize', 10, 'LineWidth', 1.5);
xlim(ppC(1) + [-hwC hwC]); ylim(ppC(2) + [-hwC hwC]); xlabel('u (px)'); ylabel('v (px)');
title('Left camera image: detected limbus points and fitted circle');
legend({'limbus points','fitted limbus','est. entry','true entry'}, ...
    'Location', 'southoutside', 'NumColumns', 2);
subplot(1,2,2);
histogram(eReg*um, 25, 'FaceColor', [0.2 0.5 0.8]); hold on;
xline(P.cam.req*um, 'r--', 'requirement', 'LineWidth', 1.5); xline(prctile95(eReg)*um, 'k:', '95th pct');
grid on; xlabel('entry-point error (\mum)'); ylabel('trials');
title(sprintf('Registration accuracy, Monte Carlo (n = %d)', nReg));
saveFig(fig, figDir, 'F10_camera_registration');

%% ======================================================================
%  S3. FORWARD / INVERSE KINEMATICS  (Modules 4-5)
%  ======================================================================
fprintf('--- S3. DH model, FK and IK ---\n');
fprintf('Standard DH table:  i | theta_i      | d_i     | a_i | alpha_i\n');
fprintf('                    1 | q1           | 0       | 0   | -90 deg\n');
fprintf('                    2 | q2           | 0       | 0   | +90 deg\n');
fprintf('                    3 | 0            | q3      | 0   |   0\n');
qT = [0 pi/2 6e-3; 0.3 1.2 4e-3; -0.5 2.0 10e-3; 0.6 0.9 -5e-3; -0.2 1.7 0.5e-3].';
fprintf('%-34s %-28s %s\n','q = [q1 q2 (deg), q3 (mm)]','tip (DH) in {0} (mm)','|DH - geometric| (m)');
errFK = zeros(1,size(qT,2));
for k = 1:size(qT,2)
    T = fp_fk(qT(:,k), P); pg = fp_fk_geom(qT(:,k));
    errFK(k) = norm(T(1:3,4) - pg);
    fprintf('[%6.1f %6.1f, %5.1f]                 [%6.3f %6.3f %6.3f]      %.1e\n', ...
        rad2deg(qT(1:2,k)), qT(3,k)*mm, T(1:3,4)*mm, errFK(k));
end
% IK candidates for site A
[qIK, ikA] = fp_ik(to0(tgtE(:,1)), P, P.q_home);
fprintf('\nIK candidates for Site A:\n  cand |   q1 (deg)  q2 (deg)  q3 (mm) | within limits | ||FK(q)-p|| (m)\n');
for k = 1:size(ikA.cand,2)
    fprintf('  %4d | %9.2f %9.2f %8.3f | %13d | %.1e\n', k, rad2deg(ikA.cand(1:2,k)), ikA.cand(3,k)*mm, ...
        ikA.within_limits(k), ikA.err(k));
end
fprintf('  selected candidate %d -> q = [%.2f deg, %.2f deg, %.3f mm]\n', ikA.selected, ...
    rad2deg(qIK(1:2)), qIK(3)*mm);
% Monte-Carlo FK(IK(p)) verification
nMC = 5000; eIK = zeros(1,nMC);
for k = 1:nMC
    qr = P.qmin + (P.qmax - P.qmin).*rand(3,1); qr(3) = 0.5e-3 + rand*15e-3;
    p = fp_fk_geom(qr); qs = fp_ik(p, P, qr);
    T = fp_fk(qs, P); eIK(k) = norm(T(1:3,4) - p);
end
fprintf('FK(IK(p)) - p over %d random targets: max %.2e m, mean %.2e m\n', nMC, max(eIK), mean(eIK));
fprintf('Deep central target -> joint-space reachable, but clearance check = %d (rejected by safety layer)\n\n', ...
    fp_lens_check(fp_ik(to0(tgtE(:,3)),P,P.q_home), P));
SUM(end+1,:) = {'Target positioning accuracy |FK(IK(P))-P|', '<= 0.1 mm', sprintf('%.1e mm', max(eIK)*mm), max(eIK) < 1e-4};

% ---- 4-DOF positioning arm (macro stage): closed-form IK, 4 branches ----
fprintf('4-DOF positioning arm (base yaw + 3-link planar arm), standard DH (m):\n');
fprintf('   i | theta | d     | a     | alpha\n');
for i = 1:4
    fprintf('   %d | q%d    | %5.3f | %5.3f | %+4.0f deg\n', i, i, P.arm.DH(i,2), P.arm.DH(i,3), rad2deg(P.arm.DH(i,4)));
end
TFdock = P.T_E0/P.arm.T_40;                         % flange pose that puts the RCM on the entry point
[qDockArm, Ida] = fp_arm_ik(TFdock, P, zeros(4,1));
fprintf('IK of the docking pose (4 candidates = front/back x elbow up/down):\n');
fprintf('  cand |   q1      q2      q3      q4  (deg) | valid\n');
for c = 1:4
    fprintf('  %4d | %7.1f %7.1f %7.1f %7.1f       | %d\n', c, rad2deg(Ida.Q(:,c)), Ida.ok(c));
end
fprintf('  selected: %s, %s -> residual %.1e m / %.1e rad\n', Ida.branch{:}, Ida.err);
rngSave = rng; rng(P.seed + 50);                    % own stream: later sections stay reproducible
nA = 3000; eA = zeros(1,nA);
for k = 1:nA
    qr = P.arm.qmin + (P.arm.qmax - P.arm.qmin).*rand(4,1);
    Tr = fp_arm_fk(qr, P); qs = fp_arm_ik(Tr, P, qr);
    eA(k) = norm(fp_arm_fk(qs, P) - Tr, 'fro');
end
rng(rngSave);
fprintf('Arm FK(IK(T)) over %d random poses: max ||T - T_d||_F = %.1e\n\n', nA, max(eA));
SUM(end+1,:) = {'Positioning-arm IK ||FK(IK(T))-T||', '<= 1e-6', sprintf('%.1e', max(eA)), max(eA) < 1e-6};
fig = figure('Color','w','Position',[50 50 1100 520]);
hold on; axis equal; grid on; box on;
cols = [0 0.45 0.74; 0.85 0.33 0.1]; t_ = linspace(0, 2*pi, 80);
hL = gobjects(1,2); nm = {'elbow up (selected)', 'elbow down (self-collision)'};
for c = [2 1]                                        % the two "front" solutions; the "back" ones turn the base by 180 deg
    [Tq, Ta] = fp_arm_fk(Ida.Q(:,c), P);
    O = [P.arm.T_Eb(1:3,4) P.arm.T_Eb(1:3,4) + [0;0;P.arm.d1] squeeze(Ta(1:3,4,2:4))]*mm;
    hL(c) = plot(O(1,:), O(3,:), '-o', 'Color', cols(c,:), 'LineWidth', 3, 'MarkerFaceColor', cols(c,:), 'MarkerSize', 7);
    pR = (Tq*P.arm.T_40)*[0;0;0;1];
    plot([O(1,end) pR(1)*mm], [O(3,end) pR(3)*mm], '--', 'Color', cols(c,:), 'LineWidth', 1.5);
end
fill(-32 + 75*cos(t_ + 0), -54 + 60*sin(t_), [0.95 0.82 0.72], 'EdgeColor', [0.7 0.55 0.45], 'FaceAlpha', 0.6);   % head section
plot(P.eye.R*cos(t_)*mm, P.eye.R*sin(t_)*mm, 'k-', 'LineWidth', 1.2);
plot(P.eye.pE(1)*mm, P.eye.pE(3)*mm, 'rp', 'MarkerFaceColor', 'r', 'MarkerSize', 13);
plot([-345 175], [-122 -122], 'k-', 'LineWidth', 3);                       % table top
text(P.arm.T_Eb(1,4)*mm + 12, 40, 'column (q_1 axis)', 'FontSize', 9);
text(135, 355, 'a_2', 'FontSize', 11, 'FontWeight', 'bold'); text(-150, 300, 'a_3', 'FontSize', 11, 'FontWeight', 'bold');
text(-345, 235, 'q_4: wrist', 'FontSize', 9);
xlabel('x_E (mm)'); ylabel('z_E (mm)'); xlim([-360 330]); ylim([-160 460]);
legend([hL, findobj(gca, 'Marker', 'p')], [nm, {'RCM target'}], 'Location', 'northeast');
title({'4-DOF arm: IK solutions of the docking pose in the arm plane', ...
    'dashed = tool unit from the flange to the RCM point;  "back" solutions (q_1 + 180^o) flip the unit and are rejected'});
saveFig(fig, figDir, 'F11_arm_ik');

%% ======================================================================
%  S4. JACOBIAN, SINGULARITIES, STATICS AND DYNAMICS  (Module 6)
%  ======================================================================
fprintf('--- S4. Jacobian, RCM velocity, singularities, statics, dynamics ---\n');
fprintf('%-30s %-16s %-18s %-16s\n','q','||J - J_num||','max lateral v_RCM','manipulability');
qJ = [qA qB [0.1;1.4;3e-3] [0.2;1.0;0.2e-3]];
for k = 1:size(qJ,2)
    [J, Jr] = fp_jacobian(qJ(:,k), P);
    Jn = zeros(3,3); hh = 1e-7;
    for i = 1:3
        e = zeros(3,1); e(i) = hh;
        Jn(:,i) = (fp_fk_geom(qJ(:,k)+e) - fp_fk_geom(qJ(:,k)-e))/(2*hh);
    end
    [~, u] = fp_fk_geom(qJ(:,k));
    lat = 0;
    for trial = 1:50                                    % random joint rates
        v = Jr*randn(3,1); lat = max(lat, norm(v - (v.'*u)*u));
    end
    Jv = J(1:3,:);
    fprintf('[%5.1f %5.1f %5.2f mm]         %-16.1e %-18.1e %-16.3e\n', rad2deg(qJ(1:2,k)), qJ(3,k)*mm, ...
        norm(J(1:3,:)-Jn), lat, abs(det(Jv)));
end
fprintf('det(Jv) = q3^2 sin(q2)  -> singular at q3 = 0 (tip at RCM) and q2 = 0/pi (needle along z0)\n');
% manipulability map
q2g = linspace(P.qmin(2), P.qmax(2), 60); q4g = linspace(0, P.qmax(3), 60);
W = zeros(numel(q4g), numel(q2g)); Ck = W;
for i = 1:numel(q4g)
    for j = 1:numel(q2g)
        J = fp_jacobian([0; q2g(j); q4g(i)], P);
        Jv = J(1:3,:);
        W(i,j) = abs(det(Jv)); Ck(i,j) = cond(Jv*diag([1 1 1e-2]));  % scale q3 to cm-equivalent
    end
end
% statics: tissue forces at the tip  tau = Jv' F
Jst = fp_jacobian(qA, P); [~, uA] = fp_fk_geom(qA);
Fax = -P.needle.F_punct*uA;                 % axial puncture reaction on the tip
nperp = cross(uA,[0;0;1]); nperp = nperp/norm(nperp);
Flat = 0.05*nperp;                          % 50 mN lateral tissue force
tau_ax = Jst(1:3,:).'*Fax; tau_lat = Jst(1:3,:).'*Flat;
fprintf('Statics at Site A: axial 0.20 N  -> tau = [%.2e %.2e N m, %.3f N]\n', tau_ax);
fprintf('                   lateral 50 mN -> tau = [%.2e %.2e N m, %.3f N]\n', tau_lat);
fprintf('  -> axial tissue loads go only to the insertion joint; lateral loads to yaw/pitch (lever arm q3)\n');
[~, ~, GA] = fp_dynamics(qA, zeros(3,1), P);
fprintf('Gravity load at Site A (static holding effort): G = [%.3f %.3f N m, %.3f N]\n', GA);

% dynamics verification
qv = [0.2; 1.4; 5e-3]; qdv = [0.3; -0.2; 2e-3];
[M, c, G, ~] = fp_dynamics(qv, qdv, P);
hh = 1e-6; Gn = zeros(3,1);
for i = 1:3
    e = zeros(3,1); e(i) = hh;
    [~,~,~,Vp] = fp_dynamics(qv+e, qdv, P); [~,~,~,Vm] = fp_dynamics(qv-e, qdv, P);
    Gn(i) = (Vp - Vm)/(2*hh);
end
% Mdot = dM/dt along qd (central difference)
Mp = fp_dynamics(qv + hh*qdv, qdv, P); Mm = fp_dynamics(qv - hh*qdv, qdv, P);
Mdot = (Mp - Mm)/(2*hh);
fprintf('Dynamics checks: min eig(M) = %.3e (>0, SPD), ||M-M''|| = %.1e\n', min(eig(M)), norm(M-M.'));
fprintf('                 ||G - dV/dq|| = %.2e,  qd''(Mdot - 2C)qd = %.2e (passivity)\n', ...
    norm(G-Gn), qdv.'*Mdot*qdv - 2*qdv.'*c);
fprintf('Reflected rotor inertia J_m N^2 = [%.2e %.2e kg m^2, %.2f kg] vs link inertia diag(M) = [%.2e %.2e, %.3f]\n', ...
    P.act.Jm.*P.act.Nr.^2, diag(M) - (P.act.Jm.*P.act.Nr.^2).');
% energy conservation in free motion (no friction, no input)
Pf = P; Pf.Bv = zeros(3);
dtE = 1e-3; tE = 0:dtE:1.0; x = [qv; qdv]; Etot = zeros(size(tE));
for k = 1:numel(tE)
    [Mk,~,~,Vk] = fp_dynamics(x(1:3), x(4:6), Pf);
    Etot(k) = 0.5*x(4:6).'*Mk*x(4:6) + Vk;
    f = @(xx) freePlant(xx, Pf);
    k1 = f(x); k2 = f(x+dtE/2*k1); k3 = f(x+dtE/2*k2); k4 = f(x+dtE*k3);
    x = x + dtE/6*(k1+2*k2+2*k3+k4);
end
relE = max(abs(Etot - Etot(1)))/abs(Etot(1));
fprintf('Energy conservation (1 s free motion, RK4): max relative drift = %.2e\n\n', relE);
SUM(end+1,:) = {'Dynamics energy drift (free motion)', '< 1e-6', sprintf('%.1e', relE), relE < 1e-6};

fig = figure('Color','w','Position',[50 50 1100 420]);
subplot(1,2,1);
contourf(rad2deg(q2g), q4g*mm, W*1e6, 20, 'LineColor','none'); colorbar;      % m^2 -> mm^2
xlabel('q_2 (deg)'); ylabel('q_3 (mm)'); title('Manipulability |det J_v| (mm^2) at q_1 = 0');
hold on; plot(rad2deg([qA(2) qA(2)]), [0 P.qmax(3)*mm], 'w--'); yline(0.5,'r-','LineWidth',1.5);
text(55, 1.5, 'q_3 \rightarrow 0: singular (tip at RCM)', 'Color','w');
subplot(1,2,2);
plot(tE, 1e10*(Etot - Etot(1))/abs(Etot(1)), 'LineWidth', 1.5); grid on;
xlabel('t (s)'); ylabel('(E(t) - E(0)) / |E(0)|   (x 10^{-10})');
title('Dynamics verification: energy conservation, free motion');
saveFig(fig, figDir, 'F02_jacobian_dynamics');

%% ======================================================================
%  S5. SENSING: ENCODER RESOLUTION REQUIREMENT  (Module 3)
%  ======================================================================
fprintf('--- S5. Encoder resolution vs tip accuracy ---\n');
res = 2*pi./(P.sens.enc_cpr*P.act.Nr(:));
lever = [P.qmax(3); P.qmax(3); 1];
tipRes = res.*lever;
names = {'J1 yaw','J2 pitch','J3 insertion'};
fprintf('%-14s %-22s %-24s\n','Joint','joint resolution','worst-case tip resolution');
for j = 1:3
    if j < 3, fprintf('%-14s %-22s %.3f um\n', names{j}, sprintf('%.2e rad (%.4f deg)', res(j), rad2deg(res(j))), tipRes(j)*um);
    else,     fprintf('%-14s %-22s %.3f um\n', names{j}, sprintf('%.3f um', res(j)*um), tipRes(j)*um); end
end
fprintf('Requirement: tip quantisation <= 10 um (10x below the 0.1 mm accuracy target) -> %s\n\n', ...
    pf(max(tipRes) <= 10e-6));
SUM(end+1,:) = {'Encoder-limited tip resolution', '<= 10 um', sprintf('%.2f um', max(tipRes)*um), max(tipRes) <= 10e-6};

%% ======================================================================
%  S6. TREMOR FILTERING, MOTION SCALING AND TELEOPERATION  (Modules 3, 7)
%  ======================================================================
fprintf('--- S6. Hand-tremor filtering and 1:10 motion scaling ---\n');
fs = 1/P.dt; tT = 0:P.dt:6; nT = numel(tT);
mj = @(s) 10*s.^3 - 15*s.^4 + 6*s.^5;
handInt = [10e-3*mj(min(tT/5,1)); 12e-3*sin(2*pi*tT/6); -8e-3*mj(min(tT/5,1))];   % intended hand motion (m)
fr = [6.0 8.5 11.0];  trem = zeros(3,nT);
for ax = 1:3
    for f0 = fr
        trem(ax,:) = trem(ax,:) + (P.trem.amp/numel(fr))*sin(2*pi*f0*tT + 2*pi*rand);
    end
end
handMeas = handInt + trem + 1e-6*randn(3,nT);                                 % 1 um sensor noise
orders = [1 2 4]; att = zeros(size(orders)); lagE = att;
skip = tT > 1;                                                                % ignore start-up
for i = 1:numel(orders)
    sos = fp_butter(orders(i), P.trem.fc, fs);
    trF = zeros(3,nT); inF = zeros(3,nT);
    for ax = 1:3
        trF(ax,:) = fp_sosfilt(sos, trem(ax,:));
        inF(ax,:) = fp_sosfilt(sos, handInt(ax,:));
    end
    att(i)  = 1 - rms(vecnorm(trF(:,skip)))/rms(vecnorm(trem(:,skip)));
    lagE(i) = P.trem.scale*rms(vecnorm(inF - handInt));
end
fprintf('Butterworth LPF, fc = %.1f Hz:  order | tremor attenuation | lag error at tip (RMS)\n', P.trem.fc);
for i = 1:numel(orders)
    fprintf('                                  %d   |      %6.2f %%     |   %.2f um\n', orders(i), 100*att(i), lagE(i)*um);
end
ordSel = 2; sos = fp_butter(ordSel, P.trem.fc, fs);
fprintf('Selected: order %d (attenuation > 90%% with the smallest lag of the passing designs)\n', ordSel);
% teleoperation command chain: scale -> filter -> virtual fixture -> IK
q0T = [deg2rad(-13); pi/2; 6e-3];  p0T = fp_fk_geom(q0T);   % start near the keep-out zone
cmdF = zeros(3,nT);
for ax = 1:3, cmdF(ax,:) = fp_sosfilt(sos, handMeas(ax,:)); end
pCmd = p0T + P.trem.scale*cmdF;  pInt = p0T + P.trem.scale*handInt;
% Virtual fixture + Cartesian speed limit (safety): the command moves toward
% the operator target by at most vTip*dt per sample; a step that would leave
% the safe region (lens / keep-out / retina / hard stop) is refused (hold).
vTip = 2e-3;                                       % m/s, max tip speed in teleoperation
QdT = zeros(3,nT); qprev = q0T; clipped = false(1,nT); pS = p0T; pSafe = zeros(3,nT);
for k = 1:nT
    stp = pCmd(:,k) - pS; ns = norm(stp);
    if ns > vTip*P.dt, stp = stp/ns*vTip*P.dt; clipped(k) = true; end
    qk = fp_ik(pS + stp, P, qprev);
    if isempty(qk) || ~fp_lens_check(qk, P) || qk(3) > P.depth_soft
        % sliding virtual fixture: remove the step component that points
        % into the forbidden region (gradient of the clearance margin)
        clipped(k) = true;
        n = marginGradient(pS, qprev, P);
        stp2 = stp - min(0, stp.'*n)*n;
        qk = fp_ik(pS + stp2, P, qprev);
        if ~isempty(qk) && fp_lens_check(qk, P) && qk(3) <= P.depth_soft
            pS = pS + stp2;
        else
            qk = fp_ik(pS, P, qprev);              % no safe direction: hold
        end
    else
        pS = pS + stp;
    end
    pSafe(:,k) = pS; QdT(:,k) = qk; qprev = qk;
end
nClip = sum(clipped);
VdT = gradient(QdT, P.dt); AdT = gradient(VdT, P.dt);
ctrlT = struct('type','CTC','wn',2*pi*8,'zeta',1,'model_scale',1,'needle',true);
ST = fp_simulate(tT, QdT, VdT, AdT, P, ctrlT);
pAct = zeros(3,nT); latV = zeros(1,nT);
for k = 1:nT
    pAct(:,k) = fp_fk_geom(ST.q(:,k));
    [~, Jr] = fp_jacobian(ST.q(:,k), P); v = Jr*ST.qd(:,k); [~,u] = fp_fk_geom(ST.q(:,k));
    latV(k) = norm(v - (v.'*u)*u);
end
eRaw   = P.trem.scale*vecnorm(trem);                 % tremor reaching the tip without filter
eTele  = vecnorm(pAct - pInt);
fprintf('Tip tremor RMS: hand %.1f um -> scaled only %.2f um -> after filter %.3f um\n', ...
    rms(vecnorm(trem))*um, rms(eRaw)*um, P.trem.scale*(1-att(orders==ordSel))*rms(vecnorm(trem))*um);
fr_ = skip & ~clipped;
fprintf('Teleoperation tip error vs intended path (fixture inactive): RMS %.2f um, max %.2f um\n', rms(eTele(fr_))*um, max(eTele(fr_))*um);
fprintf('Robot tracking error vs safe command: RMS %.2f um, max %.2f um\n', rms(vecnorm(pAct-pSafe))*um, max(vecnorm(pAct-pSafe))*um);
fprintf('Virtual-fixture interventions: %d of %d samples;  max lateral velocity at RCM = %.1e m/s\n\n', nClip, nT, max(latV));
SUM(end+1,:) = {'Tremor attenuation (RMS)', '> 90 %', sprintf('%.1f %%', 100*att(orders==ordSel)), att(orders==ordSel) > 0.9};

fig = figure('Color','w','Position',[50 50 1150 750]);
subplot(2,2,1);
for i = 1:numel(orders)
    [~, Hf, ff] = fp_butter(orders(i), P.trem.fc, fs);
    semilogx(ff, 20*log10(abs(Hf)), 'LineWidth', 1.5); hold on;
end
xline(P.trem.fband(1),'k--'); xline(P.trem.fband(2),'k--'); yline(-20,'r:');
grid on; xlim([0.1 100]); ylim([-100 5]);
xlabel('f (Hz)'); ylabel('|H| (dB)'); title('Tremor filter: magnitude response');
legend({'order 1','order 2','order 4','tremor band','','-20 dB (90 %)'},'Location','southwest');
subplot(2,2,2);
plot(tT, handMeas(2,:)*mm, 'Color',[0.7 0.7 0.7]); hold on;
plot(tT, handInt(2,:)*mm, 'k--', 'LineWidth', 1.2); plot(tT, cmdF(2,:)*mm, 'b', 'LineWidth', 1.2);
grid on; xlim([2 3]); xlabel('t (s)'); ylabel('hand y (mm)');
title('Hand input (zoom 2-3 s)'); legend({'measured','intended','filtered'});
subplot(2,2,3);
plot(tT, (pAct - p0T)*um, 'LineWidth', 1.2); hold on; plot(tT, (pInt - p0T)*um, 'k--');
yl = ylim; area(tT, yl(1) + clipped*(yl(2)-yl(1)), yl(1), 'FaceColor', [1 0.8 0.8], 'EdgeColor', 'none', 'FaceAlpha', 0.4); ylim(yl);
grid on; xlabel('t (s)'); ylabel('tip displacement in \{0\} (\mum)');
title('Teleoperated tip motion (scaled 1:10)'); legend({'x_0','y_0','z_0','intended','','','limit active'});
subplot(2,2,4);
plot(tT, eRaw*um, 'Color',[0.85 0.4 0.3]); hold on; plot(tT, vecnorm(pAct-pSafe)*um, 'b', 'LineWidth', 1.2);
grid on; xlabel('t (s)'); ylabel('tip error (\mum)');
title('Tremor at tip vs robot tracking error'); legend({'tremor, no filter','robot error (filtered)'});
saveFig(fig, figDir, 'F03_tremor_teleop');

%% ======================================================================
%  S7. MOTION PLANNING INSIDE THE EYE: C-SPACE, A*, RRT  (Module 8)
%  ======================================================================
fprintf('--- S7. Re-orientation planning (Site A -> Site B without a second puncture) ---\n');
dPlan = qA(3);
q1g = -40:1:40; q2g = 50:1:130;                       % deg grid
occ = false(numel(q2g), numel(q1g));
for i = 1:numel(q2g)
    for j = 1:numel(q1g)
        occ(i,j) = ~fp_lens_check([deg2rad(q1g(j)); deg2rad(q2g(i)); dPlan], P);
    end
end
isFree = @(qq) fp_lens_check([deg2rad(qq(1)); deg2rad(qq(2)); dPlan], P);
qs = rad2deg(qA(1:2)); qg = rad2deg(qB(1:2));
% straight joint-space line: does it collide?
nl = 200; directOK = true;
for s = linspace(0,1,nl), if ~isFree(qs + s*(qg-qs)), directOK = false; break; end, end
fprintf('Straight joint interpolation A->B collision-free: %d  (planning is required)\n', directOK);
sCell = [find(q2g==qs(2)) find(q1g==qs(1))]; gCell = [find(q2g==qg(2)) find(q1g==qg(1))];
tic; [pA, infoA] = fp_astar(occ, sCell, gCell); tA = toc;
pathA = [q1g(pA(:,2)); q2g(pA(:,1))];
rng(P.seed+1);
tic; [pathR, infoR] = fp_rrt(isFree, qs, qg, [-40;50], [40;130], ...
    struct('step',4,'goalBias',0.1,'maxIter',4000,'edgeRes',0.5,'goalTol',4)); tR = toc;
rng(P.seed+2);
pathAs = fp_shortcut(isFree, pathA, 300, 0.5);
pathRs = fp_shortcut(isFree, pathR, 300, 0.5);
plen = @(pp) sum(vecnorm(diff(pp,1,2)));
minClr = @(pp) pathClearance(pp, dPlan, P);
fprintf('%-18s %-12s %-14s %-12s %-12s %-14s\n','Planner','success','length (deg)','waypoints','time (s)','min clr (mm)');
fprintf('%-18s %-12d %-14.2f %-12d %-12.3f %-14.3f\n','A* (grid 1 deg)',~isempty(pA),plen(pathA),size(pathA,2),tA,minClr(pathA)*mm);
fprintf('%-18s %-12d %-14.2f %-12d %-12.3f %-14.3f\n','A* + shortcut',1,plen(pathAs),size(pathAs,2),NaN,minClr(pathAs)*mm);
fprintf('%-18s %-12d %-14.2f %-12d %-12.3f %-14.3f\n','RRT',infoR.success,plen(pathR),size(pathR,2),tR,minClr(pathR)*mm);
fprintf('%-18s %-12d %-14.2f %-12d %-12.3f %-14.3f\n','RRT + shortcut',infoR.success,plen(pathRs),size(pathRs,2),NaN,minClr(pathRs)*mm);
fprintf('A* expanded %d cells; RRT used %d iterations / %d nodes\n', infoA.expanded, infoR.iter, size(infoR.nodes,2));
% time-parameterise the smoothed A* path (pchip over chord length + quintic timing)
wp = pathAs; sC = [0 cumsum(vecnorm(diff(wp,1,2)))]; Lp = sC(end);
Tp = ceil(1.875*deg2rad(Lp)/(0.5*P.vmax(1))*10)/10;    % peak rate = 50 % of v_max
tP = 0:P.dt:Tp; sP = Lp*mj(tP/Tp);
q12 = [pchip(sC, wp(1,:), sP); pchip(sC, wp(2,:), sP)];
denseOK = true; clrDense = inf;
for k = 1:numel(tP)
    [okk, cl] = fp_lens_check([deg2rad(q12(:,k)); dPlan], P);
    denseOK = denseOK && okk; clrDense = min(clrDense, min(cl([1 3])));
end
fprintf('Executable trajectory: duration %.1f s, dense collision check (%d samples) = %s, min clearance %.3f mm\n\n', ...
    Tp, numel(tP), pf(denseOK), clrDense*mm);
SUM(end+1,:) = {'Planned re-orientation collision-free', 'yes, clr >= 0.2 mm', sprintf('%s, %.2f mm', pf(denseOK), clrDense*mm), denseOK};

fig = figure('Color','w','Position',[50 50 1150 470]);
subplot(1,2,1);
imagesc(q1g, q2g, occ); set(gca,'YDir','normal'); colormap(gca, [1 1 1; 0.55 0.6 0.7]); hold on;
hP = gobjects(1,5);
hP(1) = plot(pathA(1,:), pathA(2,:), 'b-', 'LineWidth', 1);
hP(2) = plot(pathAs(1,:), pathAs(2,:), 'b-', 'LineWidth', 2.5);
if infoR.success
    E = infoR.nodes; pr = infoR.parent;
    for k = 2:size(E,2), plot([E(1,k) E(1,pr(k))], [E(2,k) E(2,pr(k))], '-', 'Color', [0.9 0.7 0.7]); end
    hP(3) = plot(pathR(1,:), pathR(2,:), 'r-', 'LineWidth', 1); hP(4) = plot(pathRs(1,:), pathRs(2,:), 'r--', 'LineWidth', 2);
end
hP(5) = plot([qs(1) qg(1)], [qs(2) qg(2)], 'k:', 'LineWidth', 1.5);
plot(qs(1), qs(2), 'go', 'MarkerFaceColor','g', 'MarkerSize', 9); plot(qg(1), qg(2), 'ms', 'MarkerFaceColor','m', 'MarkerSize', 9);
xlabel('q_1 yaw (deg)'); ylabel('q_2 pitch (deg)');
title(sprintf('C-space at q_3 = %.0f mm (grey = collision)', dPlan*mm));
nmP = {'A* (grid)', 'A* + shortcut', 'RRT', 'RRT + shortcut', 'direct (collides)'};
legend(hP(isgraphics(hP)), nmP(isgraphics(hP)), 'Location', 'north', 'NumColumns', 2);
text(-35, 62, 'lens + 1 mm margin', 'FontSize', 9); text(-6, 99, 'keep-out', 'FontSize', 9);
text(-39, 127, 'retina margin', 'FontSize', 9);
subplot(1,2,2); hold on; axis equal; grid on; box on;
plot(P.eye.R*cos(tt)*mm, P.eye.R*sin(tt)*mm, 'k-');
% project onto the plane spanned by x0 and y0 through the RCM (view along z0)
for k = round(linspace(1, numel(tP), 12))
    pt = fp_fk_geom([deg2rad(q12(:,k)); dPlan]);
    plot([0 pt(2)]*mm, [0 pt(1)]*mm, '-', 'Color', [0.3 0.3 0.9 0.6], 'LineWidth', 1.5);
end
cK0 = to0(P.eye.keepout_c);
rectangle('Position', [(cK0(2)-P.eye.keepout_r)*mm (cK0(1)-P.eye.keepout_r)*mm 2*P.eye.keepout_r*mm 2*P.eye.keepout_r*mm], ...
    'Curvature', 1, 'FaceColor', [0.95 0.5 0.4]);
xlim([-5 5]); ylim([-1 9]); xlabel('y_0 (mm)'); ylabel('x_0 (mm) - insertion direction');
title('Needle sweep about the RCM (view along z_0)'); plot(0,0,'rp','MarkerFaceColor','r','MarkerSize',12);
saveFig(fig, figDir, 'F04_motion_planning');

%% ======================================================================
%  S7b. ARM TRANSPORT PLANNING AROUND THE PATIENT: C-SPACE + RRT-CONNECT  (Modules 8, 9)
%  ======================================================================
fprintf('--- S7b. Collision-free transport of the RCM unit by the 4-DOF arm (no hand-picked via points) ---\n');
CMa = fp_arm_collision_model(P); ENVa = fp_arm_env(P);
fprintf('Robot model from the 3-D geometry: %s sampled points (voxel %s mm, sampling error <= %.1f mm), %d mesh vertices for verification\n', ...
    mat2str([cellfun(@(c) size(c,2), CMa.link) size(CMa.unit,2)]), mat2str(CMa.voxel), CMa.sampErr, ...
    sum(cellfun(@(c) size(c,2), [CMa.linkFull {CMa.unitFull}])));
% the obstacle distance must never over-estimate: check against the true distance to the head mesh
[~, Vhd] = fp_head_mesh(); rngS = rng; rng(P.seed + 70);
Qc = [-150 + 330*rand(1,400); -200 + 330*rand(1,400); -60 + 280*rand(1,400)];
[dS, idS] = ENVa.sdf(Qc); dT = zeros(1,400);
for i = 1:400, dT(i) = sqrt(min(sum((Vhd - Qc(:,i).').^2, 2))); end
rng(rngS);
okH = idS == 1 & dS > 0;
fprintf('Patient distance field vs true mesh distance (%d points): over-estimate max %.2f mm (must be <= ~0), mean conservatism %.1f mm\n', ...
    nnz(okH), max(dS(okH) - dT(okH)), mean(dT(okH) - dS(okH)));
TFpre = P.T_E0/P.arm.T_40; TFpre(1:3,4) = TFpre(1:3,4) - 0.040*P.T_E0(1:3,1);   % 40 mm before the entry, on the insertion axis
qPreA = fp_arm_ik(TFpre, P, qDockArm);
qParkA = fp_arm_ik([0.34; -0.36; 0.05; -pi/2], P, deg2rad([60; 80; -100; -70]));   % parked low beside the patient's chest
verdict = {'free', 'COLLISION -> planning needed'};
sL = linspace(0, 1, 150); dLine = zeros(size(sL));
for k = 1:numel(sL), dLine(k) = fp_arm_clearance(qParkA + sL(k)*(qPreA - qParkA), P, CMa, ENVa); end
[dLmin, kL] = min(dLine); [~, iL] = fp_arm_clearance(qParkA + sL(kL)*(qPreA - qParkA), P, CMa, ENVa);
fprintf('Straight joint-space motion park -> pre-dock: min clearance %.1f mm (%s vs %s) -> %s\n', dLmin, ...
    CMa.names{iL.body}, ENVa.names{iL.obst}, verdict{1 + (dLmin <= 0)});
[sPark, sPi] = fp_arm_selfclear(qParkA, P, CMa, true); [sPre, sQi] = fp_arm_selfclear(qPreA, P, CMa, true);
fprintf('Self-collision clearance (arm capsules + unit vertices): park %.1f mm (%s), pre-dock %.1f mm (%s)\n', ...
    sPark, sPi.pair, sPre, sQi.pair);
[sUp, sUi] = fp_arm_selfclear(Ida.Q(:,1), P, CMa, true); [sDn, sDi] = fp_arm_selfclear(Ida.Q(:,2), P, CMa, true);
fprintf('Docking-pose IK branches (S3): elbow up self-clearance %.1f mm (%s) -> valid; elbow down %.1f mm (%s) -> %s\n', ...
    sUp, sUi.pair, sDn, sDi.pair, pick2(sDn <= 0, 'REJECTED (self-collision)', 'valid'));
optA = struct('dSafe', 25, 'verifyEvery', 4);
PLa = fp_arm_plan(qParkA, qPreA, P, optA);
optA.seed = 11; PLb = fp_arm_plan(qPreA, qParkA, P, optA);
for pl = {PLa, PLb; 'park -> pre-dock', 'pre-dock -> park'}
    S_ = pl{1}.stats;
    fprintf(['%-17s RRT-Connect %3d iter (%d tries), %2d -> %2d waypoints, %5.1f s | path %.2f m in %.1f s | ' ...
        'min clearance %.1f mm (%s vs %s) | min stop margin %.1f mm | self %.1f mm (%s)\n'], pl{2}, S_.iter, 3, S_.nRaw, S_.nShort, ...
        S_.time, S_.pathLenF, pl{1}.T, S_.minClr, S_.minBody, S_.minObst, S_.minSSM, S_.minSelf, S_.minSelfPair);
end
% robustness: random collision-free start poses around the table
rngS = rng; rng(P.seed + 71); nMC = 8; okMC = false(1,nMC); tMC = nan(1,nMC); cMC = nan(1,nMC); dirHit = false(1,nMC);
k = 0;
while k < nMC
    qr = P.arm.qmin + (P.arm.qmax - P.arm.qmin).*rand(4,1);
    if fp_arm_clearance(qr, P, CMa, ENVa) < 40 || fp_arm_selfclear(qr, P, CMa) < 25, continue; end
    k = k + 1;
    dl = arrayfun(@(s) fp_arm_clearance(qr + s*(qPreA - qr), P, CMa, ENVa), linspace(0, 1, 60));
    dirHit(k) = min(dl) < optA.dSafe;
    try
        Pm = fp_arm_plan(qr, qPreA, P, struct('dSafe', 25, 'nShort', 100, 'nShort2', 0, 'verifyEvery', 8, 'tries', 1, 'seed', 100 + k));
        okMC(k) = Pm.stats.okClr && Pm.stats.okSSM && Pm.stats.okSelf; tMC(k) = Pm.stats.time; cMC(k) = Pm.stats.minClr;
    catch
    end
end
rng(rngS);
nFull = nnz(cMC >= optA.dSafe);                     % paths that keep the full margin on the mesh check
fprintf(['Random start poses: %d of %d direct joint-space motions violate the %g mm clearance; collision-free plans %d / %d, ' ...
    'of which %d / %d keep the full %g mm on the mesh check (min %.1f mm; single-try planner, model sampling error <= %.1f mm), ' ...
    'planning time %.0f-%.0f s\n\n'], nnz(dirHit), nMC, optA.dSafe, nnz(okMC), nMC, nFull, nMC, optA.dSafe, min(cMC), CMa.sampErr, min(tMC), max(tMC));
armClr = min(PLa.stats.minClr, PLb.stats.minClr); armSSM = min(PLa.stats.minSSM, PLb.stats.minSSM);
armSelf = min(PLa.stats.minSelf, PLb.stats.minSelf);
SUM(end+1,:) = {'Arm transport clearance (RRT-Connect, full mesh)', sprintf('>= %g mm', optA.dSafe), sprintf('%.1f mm', armClr), armClr >= optA.dSafe};
SUM(end+1,:) = {'Arm transport stop margin (speed & separation)', '>= 0 mm', sprintf('%.1f mm', armSSM), armSSM >= 0};
SUM(end+1,:) = {'Arm self-collision clearance (transport)', sprintf('>= %g mm', PLa.stats.dSelf), sprintf('%.1f mm', armSelf), armSelf >= PLa.stats.dSelf};
SUM(end+1,:) = {'Arm planner collision-free (random start poses)', sprintf('%d / %d', nMC, nMC), sprintf('%d / %d, min %.1f mm', nnz(okMC), nMC, min(cMC)), all(okMC)};

fig = figure('Color','w','Position',[50 50 1250 520]);
subplot(1,2,1); hold on; axis equal; grid on; box on;
[Fh_, Vh_, ~, Ch_] = fp_head_mesh();
patch('Faces', Fh_, 'Vertices', Vh_, 'FaceVertexCData', Ch_, 'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceAlpha', 0.9);
rawF = zeros(3,0);
for i = 1:size(PLa.raw,2)-1
    for s_ = linspace(0, 1, 12), TF_ = fp_arm_fk(PLa.raw(:,i) + s_*(PLa.raw(:,i+1) - PLa.raw(:,i)), P); rawF(:,end+1) = TF_(1:3,4); end %#ok<SAGROW>
end
dirF = zeros(3,numel(sL));
for k = 1:numel(sL), TF_ = fp_arm_fk(qParkA + sL(k)*(qPreA - qParkA), P); dirF(:,k) = TF_(1:3,4); end
h1 = plot3(dirF(1,:)*mm, dirF(2,:)*mm, dirF(3,:)*mm, 'k--', 'LineWidth', 1.2);
h2 = plot3(rawF(1,:)*mm, rawF(2,:)*mm, rawF(3,:)*mm, '-', 'Color', [0.9 0.3 0.25], 'LineWidth', 1);
h3 = plot3(PLa.pF(1,:)*mm, PLa.pF(2,:)*mm, PLa.pF(3,:)*mm, '-', 'Color', [0 0.45 0.8], 'LineWidth', 2.5);
h4 = plot3(PLb.pF(1,:)*mm, PLb.pF(2,:)*mm, PLb.pF(3,:)*mm, '-', 'Color', [0.1 0.65 0.3], 'LineWidth', 2);
view(55, 25); camlight; lighting gouraud; xlabel('x_E (mm)'); ylabel('y_E (mm)'); zlabel('z_E (mm)');
legend([h1 h2 h3 h4], {sprintf('direct (collides, %.0f mm)', dLmin), 'raw RRT', ...
    'planned', 'return'}, 'Location', 'northoutside');
title('Flange paths of the 4-DOF arm around the patient');
subplot(1,2,2);
yyaxis left; plot(PLa.t, PLa.clr, 'LineWidth', 1.5); hold on; plot(PLa.t(~isnan(PLa.clrFull)), PLa.clrFull(~isnan(PLa.clrFull)), '.');
yline(optA.dSafe, '--', 'd_{safe}'); ylabel('clearance lower bound (mm)');
yyaxis right; plot(PLa.t, max(PLa.v, PLa.vR)*mm, 'LineWidth', 1.5); plot(PLa.t, PLa.vlim*mm, ':', 'LineWidth', 1.5);
ylabel('speed (mm/s)'); xlabel('t (s)'); grid on;
legend({'clearance (model)', 'clearance (mesh)', '', 'speed', 'SSM limit'}, 'Location', 'southoutside', 'NumColumns', 2);
title({'Speed and separation monitoring along the planned path', 'v \leq v_{lim}(d):  v T_r + v^2/(2a) + margin \leq clearance'});
saveFig(fig, figDir, 'F13_arm_planning');

%% ======================================================================
%  S8. COMPLETE INJECTION MISSION: TRAJECTORIES + CONTROL  (Modules 3, 7)
%  ======================================================================
fprintf('--- S8. Mission: align -> approach -> insert -> inject A -> re-orient -> inject B -> retract ---\n');
dt = P.dt; Qs = []; Vs = []; As = []; phase = [];
optT = struct('vmax', P.vmax, 'amax', P.amax);
qAout = [qA(1:2); -10e-3];
[Qs,Vs,As,phase] = addSeg(Qs,Vs,As,phase,'trapezoid',P.q_home,qAout,optT,dt,1);   % 1 align
qAsurf = [qA(1:2); -0.5e-3];
[Qs,Vs,As,phase] = addSeg(Qs,Vs,As,phase,'trapezoid',qAout,qAsurf,optT,dt,2);      % 2 approach
[Qs,Vs,As,phase] = addSeg(Qs,Vs,As,phase,'quintic',qAsurf,qA,struct('T',6),dt,3);  % 3 insertion
[Qs,Vs,As,phase] = addHold(Qs,Vs,As,phase,qA,2.0,dt,4);                            % 4 inject A
q12r = deg2rad(q12); Q5 = [q12r; dPlan*ones(1,numel(tP))];
V5 = gradient(Q5, dt); A5 = gradient(V5, dt);
Qs = [Qs Q5(:,2:end)]; Vs = [Vs V5(:,2:end)]; As = [As A5(:,2:end)]; phase = [phase 5*ones(1,numel(tP)-1)]; % 5 re-orient
[Qs,Vs,As,phase] = addHold(Qs,Vs,As,phase,qB,2.0,dt,6);                            % 6 inject B
qBout = [qB(1:2); -10e-3];
[Qs,Vs,As,phase] = addSeg(Qs,Vs,As,phase,'quintic',qB,qBout,struct('T',7),dt,7);   % 7 retract
tM = (0:size(Qs,2)-1)*dt;
fprintf('Mission duration %.1f s, %d samples at %.0f Hz\n', tM(end), numel(tM), 1/dt);
fprintf('Peak desired rates: |q1d| %.1f deg/s, |q2d| %.1f deg/s, |q3d| %.2f mm/s (limits %.0f, %.0f deg/s, %.0f mm/s)\n', ...
    rad2deg(max(abs(Vs(1,:)))), rad2deg(max(abs(Vs(2,:)))), max(abs(Vs(3,:)))*mm, ...
    rad2deg(P.vmax(1)), rad2deg(P.vmax(2)), P.vmax(3)*mm);

ctrls = {struct('type','PID','wn',2*pi*8,'zeta',1,'model_scale',1,'needle',true), ...
         struct('type','CTC','wn',2*pi*8,'zeta',1,'model_scale',1,'needle',true), ...
         struct('type','CTC','wn',2*pi*8,'zeta',1,'model_scale',0.8,'needle',true)};
cname = {'PID','CTC exact','CTC -20 % mass'};
SM = cell(1,3); tipE = zeros(3,numel(tM));
pDes = zeros(3,numel(tM)); for k = 1:numel(tM), pDes(:,k) = fp_fk_geom(Qs(:,k)); end
inside = Qs(3,:) > 0;
fprintf('%-24s %-14s %-14s %-16s %-14s %-10s\n','Controller','tip RMS (um)','tip max (um)','tip max inside','|e_q3| max','sat %');
for c = 1:3
    SM{c} = fp_simulate(tM, Qs, Vs, As, P, ctrls{c});
    pa = zeros(3,numel(tM)); for k = 1:numel(tM), pa(:,k) = fp_fk_geom(SM{c}.q(:,k)); end
    tipE(c,:) = vecnorm(pa - pDes);
    sat = 100*mean(any(abs(SM{c}.tau_cmd) > SM{c}.tau_max + 1e-12, 1));
    fprintf('%-24s %-14.2f %-14.2f %-16.2f %-14.2f %-10.2f\n', cname{c}, rms(tipE(c,:))*um, max(tipE(c,:))*um, ...
        max(tipE(c,inside))*um, max(abs(SM{c}.e(3,:)))*um, sat);
end
SUM(end+1,:) = {'Mission tip tracking error (CTC, inside eye)', '<= 100 um', ...
    sprintf('%.1f um max', max(tipE(2,inside))*um), max(tipE(2,inside)) <= 1e-4};
SUM(end+1,:) = {'Mission tip error with -20 % model error', '<= 100 um', ...
    sprintf('%.1f um max', max(tipE(3,inside))*um), max(tipE(3,inside)) <= 1e-4};
% RCM lateral motion during the whole mission (CTC exact)
latM = zeros(1,numel(tM));
for k = 1:numel(tM)
    [~, Jr] = fp_jacobian(SM{2}.q(:,k), P); v = Jr*SM{2}.qd(:,k); [~,u] = fp_fk_geom(SM{2}.q(:,k));
    latM(k) = norm(v - (v.'*u)*u);
end
rcmDrift = sum(latM)*dt;
fprintf('Ideal mechanism: integrated lateral motion of the needle at the RCM = %.2e m (numerical zero)\n', rcmDrift);

% puncture detection from the noisy force sensor
Ftrue = SM{2}.F; Fm = Ftrue + P.sens.force_sd*randn(size(Ftrue));
Ff = filter(ones(1,5)/5, 1, Fm);                       % 5 ms moving average
runMax = 0; kDet = NaN;
for k = 6:numel(Ff)
    runMax = max(runMax, Ff(k));
    if runMax > 0.5*P.needle.F_punct && Ff(k) < 0.4*runMax && SM{2}.qd(3,k) > 0
        kDet = k; break;
    end
end
kTrue = find(SM{2}.punct, 1);
fprintf('Scleral puncture: true at t = %.3f s (q3 = %.2f mm), detected at t = %.3f s -> delay %.1f ms\n\n', ...
    tM(kTrue), SM{2}.q(3,kTrue)*mm, tM(kDet), (tM(kDet)-tM(kTrue))*1e3);
SUM(end+1,:) = {'Puncture detection delay (force sensor)', '<= 20 ms', sprintf('%.1f ms', (tM(kDet)-tM(kTrue))*1e3), (tM(kDet)-tM(kTrue)) <= 0.02};

fig = figure('Color','w','Position',[50 50 1150 780]);
lbl = {'q_1 (deg)','q_2 (deg)','q_3 (mm)'}; sc = [180/pi 180/pi mm];
for j = 1:3
    subplot(3,2,2*j-1);
    plot(tM, Qs(j,:)*sc(j), 'k--', 'LineWidth', 1.5); hold on;
    plot(tM, SM{2}.q(j,:)*sc(j), 'b', 'LineWidth', 1);
    grid on; ylabel(lbl{j}); if j == 1, title('Joint tracking (CTC exact)'); legend({'desired','actual'}); end
    phaseLines(tM, phase);
end
xlabel('t (s)');
subplot(3,2,2);
semilogy(tM, max(tipE,1e-9)*um); grid on; ylim([1e-2 1e3]); ylabel('tip error (\mum)'); title('Tip tracking error');
legend(cname, 'Location','southeast'); phaseLines(tM, phase);
subplot(3,2,4);
plot(tM, Fm, 'Color', [0.7 0.7 0.7]); hold on; plot(tM, Ftrue, 'b', 'LineWidth', 1.2);
xline(tM(kTrue), 'g-', 'puncture', 'LabelHorizontalAlignment','left'); xline(tM(kDet), 'r--', 'detected', 'LabelHorizontalAlignment','right');
grid on; ylabel('needle force (N)'); title('Needle force sensing and puncture detection');
xlim([tM(kTrue)-1.5 tM(kTrue)+1.5]); legend({'measured','true'});
subplot(3,2,6);
plot(tM, SM{2}.tau(1:2,:)*1e3, 'LineWidth', 1); hold on;
yyaxis right; plot(tM, SM{2}.tau(3,:), 'LineWidth', 1); ylabel('F_3 (N)');
yyaxis left; ylabel('\tau_{1,2} (mN m)'); grid on; xlabel('t (s)');
title('Joint efforts (CTC exact)'); legend({'\tau_1','\tau_2','F_3'});
saveFig(fig, figDir, 'F05_mission_control');

%% ======================================================================
%  S9. ACTUATOR AND END-EFFECTOR FEASIBILITY  (Module 2)
%  ======================================================================
fprintf('--- S9. Drive feasibility from the simulated mission (CTC exact) ---\n');
A = P.act; tauJ = SM{2}.tau; qdJ = SM{2}.qd;
tm = zeros(3,2); wmx = zeros(3,1);
for j = 1:2
    tm(j,:) = [max(abs(tauJ(j,:))) rms(tauJ(j,:))]/(A.N(j)*A.eta(j));
    wmx(j) = max(abs(qdJ(j,:)))*A.N(j);
end
k4 = A.lead/(2*pi*A.eta_s*A.N4);                       % N -> motor N m
tm(3,:) = [max(abs(tauJ(3,:))) rms(tauJ(3,:))]*k4;
wmx(3) = max(abs(qdJ(3,:)))*2*pi*A.N4/A.lead;
Trated = [A.Tm_rated A.Tm4_rated]; Wrated = [A.wm_rated A.wm4_rated];
fprintf('%-14s %-14s %-14s %-12s %-14s %-12s %s\n','Joint','peak Tm (mNm)','RMS Tm (mNm)','rated(mNm)','speed (rpm)','rated (rpm)','status (SF 1.5)');
okDrive = true;
for j = 1:3
    okj = A.SF*tm(j,2) <= Trated(j) && tm(j,1) <= 3*Trated(j) && wmx(j) <= Wrated(j);
    okDrive = okDrive && okj;
    fprintf('%-14s %-14.3f %-14.3f %-12.1f %-14.0f %-12.0f %s\n', names{j}, tm(j,1)*1e3, tm(j,2)*1e3, Trated(j)*1e3, ...
        wmx(j)*60/(2*pi), Wrated(j)*60/(2*pi), pf(okj));
end
% emergency retraction demand on the insertion drive
mEff = A.Jm(3)*A.Nr(3)^2 + P.link(3).m;
Fret = mEff*P.safe.a_retract + P.Bv(3,3)*P.safe.v_retract + P.needle.f_fric*P.qmax(3);
Tret = Fret*k4; wret = P.safe.v_retract*2*pi*A.N4/A.lead;
fprintf('Emergency retraction (a = %.0f m/s^2, v = %.0f mm/s): effective mass %.2f kg, F = %.1f N -> Tm = %.1f mN m (peak avail. %.1f), %.0f rpm (rated %.0f)\n', ...
    P.safe.a_retract, P.safe.v_retract*mm, mEff, Fret, Tret*1e3, 3*A.Tm4_rated*1e3, wret*60/(2*pi), A.wm4_rated*60/(2*pi));
etaBack = 2 - 1/A.eta_s;
fprintf('Lead-screw back-drive efficiency = 2 - 1/eta = %.2f (< 0 -> self-locking): needle cannot be pushed deeper when power is lost\n', etaBack);
fprintf('End effector: 30G needle (OD %.3f mm, ID %.3f mm, L %.1f mm) + 1 mL syringe driven by a micro syringe pump\n\n', ...
    P.needle.OD*mm, P.needle.ID*mm, P.needle.L*mm);
SUM(end+1,:) = {'Motor torque/speed feasibility (all joints)', 'pass with SF 1.5', pf(okDrive), okDrive};
SUM(end+1,:) = {'Insertion drive self-locking (power loss)', 'back-drive eta < 0', sprintf('%.2f', etaBack), etaBack < 0};

fig = figure('Color','w','Position',[50 50 900 380]);
bar(categorical(names, names), [tm(:,1)./Trated(:) tm(:,2)./Trated(:)]); hold on;
yline(1,'r--','rated'); yline(3,'r:','peak (3x)','LabelHorizontalAlignment','left'); grid on; ylim([0 3.4]);
ylabel('motor torque / rated torque'); legend({'peak','RMS'});
title('Actuator utilisation along the mission (CTC exact)');
saveFig(fig, figDir, 'F06_actuators');

%% ======================================================================
%  S9b. 4-DOF ARM: DYNAMICS, DRIVE SIZING, BRAKES, BRAKING CAPABILITY  (Modules 2, 6, 9)
%  ======================================================================
fprintf('--- S9b. Positioning-arm dynamics and drive sizing (arm + RCM unit %.2f kg) ---\n', P.arm.payload.m);
fprintf('%-8s %-8s %-26s %-30s\n', 'Body', 'm (kg)', 'COM in own frame (mm)', 'principal inertia (kg m^2)');
bn = {'link 1', 'link 2', 'link 3', 'link 4', 'RCM unit'}; bb = [P.arm.link, P.arm.payload];
for i = 1:5
    fprintf('%-8s %-8.2f [%6.1f %6.1f %6.1f]       [%s]\n', bn{i}, bb(i).m, bb(i).rc*mm, num2str(sort(eig(bb(i).I)).', '%.2e '));
end
qDk = qDockArm;
[Ma, ca, Ga] = fp_arm_dynamics(qDk, [0.3; -0.2; 0.4; 0.5], P, false);
fprintf('Docked pose q = [%s] deg:\n', num2str(rad2deg(qDk).', '%.1f '));
fprintf('  M (links + unit) = [%s; %s; %s; %s] kg m^2\n', num2str(Ma(1,:), '%.4f '), num2str(Ma(2,:), '%.4f '), ...
    num2str(Ma(3,:), '%.4f '), num2str(Ma(4,:), '%.4f '));
fprintf('  G = [%s] N m  (gravity holding torque while the brakes are locked)\n', num2str(Ga.', '%.3f '));
hh = 1e-6; Gn = zeros(4,1); qdv = [0.3; -0.2; 0.4; 0.5];
for i = 1:4
    e = zeros(4,1); e(i) = hh;
    [~,~,~,Vp] = fp_arm_dynamics(qDk + e, [], P, false); [~,~,~,Vm] = fp_arm_dynamics(qDk - e, [], P, false);
    Gn(i) = (Vp - Vm)/(2*hh);
end
Mp = fp_arm_dynamics(qDk + hh*qdv, [], P, false); Mm = fp_arm_dynamics(qDk - hh*qdv, [], P, false);
fprintf('  checks: min eig(M) = %.3e (SPD), ||G - dV/dq|| = %.1e, qd''(Mdot - 2C)qd = %.1e\n', min(eig(Ma)), norm(Ga - Gn), ...
    qdv.'*((Mp - Mm)/(2*hh))*qdv - 2*qdv.'*ca);
% worst static gravity torque anywhere in the joint space (brake / holding requirement)
rngS = rng; rng(P.seed + 80); nG = 20000; Gw = zeros(4,1); qGw = zeros(4,4);
for k = 1:nG
    qr = P.arm.qmin + (P.arm.qmax - P.arm.qmin).*rand(4,1);
    [~,~,Gk] = fp_arm_dynamics(qr, [], P, false);
    for j = 1:4, if abs(Gk(j)) > Gw(j), Gw(j) = abs(Gk(j)); qGw(:,j) = qr; end, end
end
rng(rngS);
fprintf('Worst static gravity torque over %d random poses: |G|max = [%s] N m\n', nG, num2str(Gw.', '%.2f '));
% joint motion along both planned transport paths (S7b) + rotor-free torques
paths = {PLa, PLb}; QQ = []; QD = []; QDD = []; TL_ = [];
for i = 1:2
    dtp = paths{i}.t(2) - paths{i}.t(1); Q_ = paths{i}.Q; Qd_ = gradient(Q_, dtp); Qdd_ = gradient(Qd_, dtp);
    for j = 1:4, Qdd_(j,:) = conv(Qdd_(j,:), ones(1,5)/5, 'same'); end
    QQ = [QQ Q_]; QD = [QD Qd_]; QDD = [QDD Qdd_]; TL_ = [TL_ paths{i}.t + (i-1)*(paths{1}.T + 1)]; %#ok<AGROW>
end
nP = size(QQ,2); tauL = zeros(4,nP);
for k = 1:nP
    [Mk, ck, Gk] = fp_arm_dynamics(QQ(:,k), QD(:,k), P, false); tauL(:,k) = Mk*QDD(:,k) + ck + Gk;
end
fprintf('Joint demand on the transport paths: peak |tau| = [%s] N m, RMS = [%s] N m, peak |qd| = [%s] rad/s\n', ...
    num2str(max(abs(tauL), [], 2).', '%.2f '), num2str(sqrt(mean(tauL.^2, 2)).', '%.2f '), num2str(max(abs(QD), [], 2).', '%.2f '));
% sizing workflow (M2): smallest motor class + ratio that meets SF*T_rms <= T_rated, T_peak <= T_pk,
% motor speed <= rated, and the worst static pose can be held by the motor at peak torque
cat_ = struct('name', {'30 W', '70 W', '150 W', '250 W'}, 'Tr', {0.05, 0.10, 0.25, 0.45}, 'Tp', {0.15, 0.30, 0.75, 1.35}, ...
    'wr', {6000*2*pi/60, 5000*2*pi/60, 4500*2*pi/60, 4000*2*pi/60}, 'Jm', {4e-6, 1.2e-5, 3.3e-5, 7.0e-5});
ratios = [100 120 160]; etaH = 0.75; SFa = P.act.SF; sel = zeros(4,2);
fprintf('%-5s %-7s %-6s %-12s %-12s %-12s %-12s %-10s %s\n', 'Joint', 'motor', 'N', 'T_rms*SF', 'T_rated', 'T_peak', 'T_pk avail', 'rpm', 'hold worst pose');
okArm = true; Tsel = zeros(4,6);
for j = 1:4
    found = false;
    for c = 1:numel(cat_)
        for N_ = ratios
            tm = (abs(tauL(j,:)) + cat_(c).Jm*N_^2*abs(QDD(j,:)))/(N_*etaH);
            trms = sqrt(mean(tm.^2)); tpk = max(tm); wmx = max(abs(QD(j,:)))*N_; thold = Gw(j)/(N_*etaH);
            if SFa*trms <= cat_(c).Tr && tpk <= cat_(c).Tp && wmx <= cat_(c).wr && thold <= cat_(c).Tp
                sel(j,:) = [c N_]; Tsel(j,:) = [SFa*trms cat_(c).Tr tpk cat_(c).Tp wmx*60/(2*pi) thold]; found = true; break
            end
        end
        if found, break; end
    end
    okArm = okArm && found;
    if found
        fprintf('J%-4d %-7s %-6d %-12.4f %-12.3f %-12.4f %-12.3f %-10.0f %.3f N m (motor side)\n', j, cat_(sel(j,1)).name, sel(j,2), Tsel(j,1:5), Tsel(j,6));
    else
        fprintf('J%-4d no drive in the catalogue meets the demand\n', j);
    end
end
% use the selected drives from here on (reflected inertia, braking)
for j = 1:4
    if sel(j,1) > 0
        P.arm.act.N(j) = sel(j,2); P.arm.act.Tm_rated(j) = cat_(sel(j,1)).Tr; P.arm.act.Tm_peak(j) = cat_(sel(j,1)).Tp;
        P.arm.act.wm_rated(j) = cat_(sel(j,1)).wr; P.arm.act.Jm(j) = cat_(sel(j,1)).Jm;
    end
end
P.arm.act.Tbrake = 2*Gw.';                                           % power-off brakes: 2 x worst gravity torque
fprintf('Power-off brakes (spring-applied): holding torque = 2 x worst gravity = [%s] N m at the joint = [%s] N m at the motor\n', ...
    num2str(P.arm.act.Tbrake, '%.1f '), num2str(P.arm.act.Tbrake./P.arm.act.N, '%.3f '));
fprintf('  brake check at the docked pose: |G| / T_brake = [%s] (< 1 -> arm stays locked during the injection)\n', ...
    num2str((abs(Ga.')./max(P.arm.act.Tbrake, eps)).*(P.arm.act.Tbrake > 0), '%.2f '));
% protective-stop braking capability along the planned paths (replaces the assumed 1 m/s^2)
tauPk = P.arm.act.Tm_peak.*P.arm.act.N.*P.arm.act.eta;              % joint torque available with peak motor torque
aCap = inf(1,nP);
for k = 1:nP
    TF_ = fp_arm_fk(QQ(:,k), P); dqk = QD(:,k);
    [~, Ta_] = fp_arm_fk(QQ(:,k), P); o_ = [P.arm.T_Eb(1:3,4) squeeze(Ta_(1:3,4,1:3))]; z_ = [P.arm.T_Eb(1:3,3) squeeze(Ta_(1:3,3,1:3))];
    Jp = zeros(3,4); for j = 1:4, Jp(:,j) = cross(z_(:,j), TF_(1:3,4) - o_(:,j)); end
    vF = norm(Jp*dqk); if vF < 0.02, continue; end                   % moving at > 20 mm/s only
    tq = dqk/vF;                                                     % joint motion per metre of flange path
    [Mk, ck, Gk] = fp_arm_dynamics(QQ(:,k), dqk, P, true);
    mi = Mk*tq; hi = ck + Gk;                                        % tau(s_dd) = hi - s_dd * mi
    bnd = (tauPk(:) - sign(mi).*hi)./abs(mi); bnd(abs(mi) < 1e-9) = inf;
    aCap(k) = min(bnd);
end
aCapMin = min(aCap);
fprintf(['Protective stop (cat. 1, motors braking at peak torque): minimum flange deceleration along the paths = %.2f m/s^2 ' ...
    '(assumed in S13: %.2f m/s^2) -> %s\n\n'], aCapMin, P.safe.arm_abrake, pf(aCapMin >= P.safe.arm_abrake));
SUM(end+1,:) = {'Arm drive torque/speed feasibility (all joints)', 'pass with SF 1.5', pf(okArm), okArm};
SUM(end+1,:) = {'Arm braking capability vs assumed deceleration', sprintf('>= %.1f m/s^2', P.safe.arm_abrake), sprintf('%.2f m/s^2', aCapMin), aCapMin >= P.safe.arm_abrake};

fig = figure('Color','w','Position',[50 50 1250 460]);
subplot(1,2,1); hold on; grid on; box on;
plot(TL_, tauL.', 'LineWidth', 1.4); xlabel('t (s)  -  park -> pre-dock, then pre-dock -> park'); ylabel('joint torque (N m)');
legend({'\tau_1 base', '\tau_2 shoulder', '\tau_3 elbow', '\tau_4 wrist'}, 'Location', 'best');
title({'4-DOF arm: joint torques along the', 'planned transport paths (with RCM unit)'});
subplot(1,2,2);
Y = [Tsel(:,1)./Tsel(:,2) Tsel(:,3)./Tsel(:,4) Tsel(:,6)./Tsel(:,4)];
bar(categorical({'J1', 'J2', 'J3', 'J4'}), Y); hold on; yline(1, 'r--', 'limit'); grid on;
ylabel('utilisation'); legend({'SF x RMS / rated', 'peak / peak rating', 'static hold / peak'}, 'Location', 'northoutside');
title({'Selected drives', strjoin(arrayfun(@(j) sprintf('J%d %s 1:%d', j, cat_(max(sel(j,1),1)).name, sel(j,2)), 1:4, 'UniformOutput', false), ', ')});
saveFig(fig, figDir, 'F14_arm_drives');

%% ======================================================================
%  S10. RCM ACCURACY UNDER MANUFACTURING ERRORS + PIVOT CALIBRATION  (M4-M6, M9)
%  ======================================================================
fprintf('--- S10. RCM drift under DH parameter errors (Monte Carlo) ---\n');
nRob = 300; nCfg = 40;
sig0 = struct('th', deg2rad(0.05), 'd', 20e-6, 'a', 20e-6, 'al', deg2rad(0.05));
[dNom, dCal] = rcmMonteCarlo(P, sig0, nRob, nCfg, 1);
fprintf('Errors (1 sigma): theta/alpha %.2f deg, a/d %.0f um, %d robots x %d poses\n', ...
    rad2deg(sig0.th), sig0.d*um, nRob, nCfg);
fprintf('RCM drift before calibration: median %.1f um, 95th pct %.1f um\n', median(dNom)*um, prctile95(dNom)*um);
fprintf('RCM drift after pivot calibration (least-squares RCM): median %.1f um, 95th pct %.1f um\n', ...
    median(dCal)*um, prctile95(dCal)*um);
scl = [0.1 0.2 0.3 0.5 0.75 1 1.5 2]; p95 = zeros(size(scl)); p95c = p95;
for i = 1:numel(scl)
    sg = struct('th', scl(i)*sig0.th, 'd', scl(i)*sig0.d, 'a', scl(i)*sig0.a, 'al', scl(i)*sig0.al);
    [dn, dc] = rcmMonteCarlo(P, sg, 150, 30, 10+i);
    p95(i) = prctile95(dn); p95c(i) = prctile95(dc);
end
tolOK = scl(find(p95c <= 10e-6, 1, 'last'));
if isempty(tolOK), tolOK = NaN; p95tol = NaN; else, p95tol = p95c(scl == tolOK); end
fprintf('95th-pct calibrated drift vs tolerance scale [%s]: [%s] um\n', num2str(scl), num2str(p95c*um, '%7.1f'));
fprintf('Largest tolerance scale with 95 %% of calibrated robots below 10 um: %.2f x (theta/alpha %.3f deg, a/d %.1f um)\n', ...
    tolOK, rad2deg(tolOK*sig0.th), tolOK*sig0.d*um);
fprintf('-> with the nominal tolerances the calibrated drift does NOT meet 10 um; the tolerance above becomes the design specification\n\n');
SUM(end+1,:) = {'RCM drift, nominal tolerances (calibrated, 95th pct)', '< 10 um', sprintf('%.1f um', prctile95(dCal)*um), prctile95(dCal) < 10e-6};
SUM(end+1,:) = {'RCM drift, specified tolerances (calibrated, 95th pct)', '< 10 um', sprintf('%.1f um @ %.2fx', p95tol*um, tolOK), p95tol < 10e-6};

fig = figure('Color','w','Position',[50 50 1100 400]);
subplot(1,2,1);
histogram(dNom*um, 30, 'FaceColor', [0.85 0.4 0.3], 'FaceAlpha', 0.6); hold on;
histogram(dCal*um, 30, 'FaceColor', [0.2 0.5 0.8], 'FaceAlpha', 0.6); xline(10,'k--','10 \mum');
xlabel('max RCM drift over the pose set (\mum)'); ylabel('robots');
legend({'nominal RCM','after pivot calibration'}); title('RCM drift from assembly errors'); grid on;
subplot(1,2,2);
plot(scl, p95*um, 'o-', scl, p95c*um, 's-', 'LineWidth', 1.5); yline(10,'k--','requirement');
grid on; xlabel('tolerance scale (x nominal \sigma)'); ylabel('95th percentile drift (\mum)');
legend({'uncalibrated','calibrated'},'Location','northwest'); title('Tolerance allocation');
saveFig(fig, figDir, 'F07_rcm_calibration');

%% ======================================================================
%  S11. FLUIDIC INJECTION AND INTRAOCULAR PRESSURE  (Modules 2, 3, 7, 9)
%  ======================================================================
fprintf('--- S11. Drug delivery (50 uL) and intraocular pressure ---\n');
% verification: closed eye, instantaneous bolus -> P = P0 exp(K V)
Pv = P; Pv.iop.C = 0; Pv.iop.Fu = 0;
Rv = fp_iop_sim('bolus', Pv, struct('Tend', 2, 'Tbolus', 0.5, 'margin', 0, 'noise', false));
Pan = P.iop.P0*exp(P.iop.K*P.iop.Vdose);
fprintf('Verification (no outflow): simulated %.3f mmHg vs analytical P0*exp(K*V) = %.3f mmHg (err %.2e)\n', ...
    Rv.P(end), Pan, abs(Rv.P(end)-Pan));
optI = struct('Tend', 1800, 'Tbolus', 3, 'Tconst', 60, 'margin', 1, 'Kp', 0.5, 'Qcap', 5, 'noise', true);
strat = {'bolus','constant','feedback'}; sname = {'bolus 3 s','constant 60 s','IOP feedback'};
RI = cell(1,3);
fprintf('%-26s %-11s %-12s %-14s %-14s %-12s %-10s %s\n','Strategy','peak IOP','t_inject','t>25mmHg (s)','recovery (min)','dP needle','F plunger','Re');
for s = 1:3
    rng(P.seed+20);
    RI{s} = fp_iop_sim(strat{s}, P, optI);
    fprintf('%-26s %-11.1f %-12.1f %-14.1f %-14.1f %-12.1f %-10.3f %.1f\n', sname{s}, RI{s}.Ppeak, RI{s}.Tinj, ...
        RI{s}.time_above_lim, RI{s}.Trecover/60, RI{s}.dPmax_needle, RI{s}.Fplunger, RI{s}.Re);
end
fprintf('(peak IOP mmHg, t_inject s, needle pressure drop mmHg, plunger force N; Re << 2000 -> laminar)\n');
SUM(end+1,:) = {'Peak IOP (robot feedback injection)', '< 25 mmHg', sprintf('%.1f mmHg', RI{3}.Ppeak), RI{3}.Ppeak < 25};
% Monte Carlo over patients (ocular rigidity and baseline IOP)
nPat = 120; Kp_ = 0.008 + 0.014*rand(1,nPat); P0p = min(max(15 + 2.5*randn(1,nPat), 10), 21);
pk = zeros(3,nPat); ti = zeros(3,nPat);
optMC = optI; optMC.Tend = 900;
for i = 1:nPat
    for s = 1:3
        o = optMC; o.K = Kp_(i); o.P0 = P0p(i);
        Ri = fp_iop_sim(strat{s}, P, o); pk(s,i) = Ri.Ppeak; ti(s,i) = Ri.Tinj;
    end
end
fprintf('Patient Monte Carlo (K in [0.008 0.022] 1/uL, P0 ~ N(15,2.5) mmHg, n = %d):\n', nPat);
for s = 1:3
    fprintf('  %-24s peak IOP: median %.1f, max %.1f mmHg | %% patients > 25 mmHg: %.0f %% | > 40 mmHg: %.0f %% | t_inj median %.0f s\n', ...
        sname{s}, median(pk(s,:)), max(pk(s,:)), 100*mean(pk(s,:) > 25), 100*mean(pk(s,:) > 40), median(ti(s,:),'omitnan'));
end
% trade-off: allowed IOP limit vs injection time
PlimS = [22 25 28 31 35 40]; TinjS = zeros(size(PlimS));
for i = 1:numel(PlimS)
    Pt = P; Pt.iop.Plim = PlimS(i);
    Rt = fp_iop_sim('feedback', Pt, struct('Tend', 900, 'margin', 1, 'Kp', 0.5, 'Qcap', 5, 'noise', false));
    TinjS(i) = Rt.Tinj;
end
fprintf('Trade-off (nominal eye): IOP limit [%s] mmHg -> injection time [%s] s\n\n', num2str(PlimS), num2str(round(TinjS)));

fig = figure('Color','w','Position',[50 50 1150 750]);
cols = [0.85 0.33 0.1; 0.47 0.67 0.19; 0 0.45 0.74];
subplot(2,2,1); hold on;
for s = 1:3, plot(RI{s}.t/60, RI{s}.P, 'Color', cols(s,:), 'LineWidth', 1.5); end
yline(P.iop.Plim,'k--','25 mmHg target'); yline(P.iop.Palarm,'r:','alarm');
grid on; xlabel('t (min)'); ylabel('IOP (mmHg)'); title('Intraocular pressure after a 50 \muL injection');
legend(sname, 'Location', 'northeast'); xlim([0 30]);
subplot(2,2,2); hold on;
hV = gobjects(1,3);                        % plot only the displayed window (no lines outside the axes)
for s = 1:3, kk = RI{s}.t <= 360; hV(s) = plot(RI{s}.t(kk), RI{s}.V(kk), '-', 'Color', cols(s,:), 'LineWidth', 1.5, 'Marker', 'none'); end
grid on; xlabel('t (s)'); ylabel('injected volume (\muL)'); title('Delivered dose');
yyaxis right; kk = RI{3}.t <= 360;
hQ = plot(RI{3}.t(kk), RI{3}.Q(kk), '-', 'Color', [0.45 0.70 0.90], 'Marker', 'none'); ylabel('feedback flow (\muL/s)'); xlim([0 360]);
legend([hV hQ], [sname, {'feedback flow (right axis)'}], 'Location', 'east');
subplot(2,2,3);
hold on;                                   % jittered scatter + median (base MATLAB)
for s = 1:3
    scatter(s + 0.15*randn(1,nPat), pk(s,:), 10, cols(s,:), 'filled', 'MarkerFaceAlpha', 0.4);
    plot(s + [-0.3 0.3], median(pk(s,:))*[1 1], 'k-', 'LineWidth', 2);
end
set(gca, 'XTick', 1:3, 'XTickLabel', {'bolus','constant','feedback'}); xlim([0.4 3.6]);
yline(25,'k--'); yline(40,'r:');
ylabel('peak IOP (mmHg)'); title(sprintf('Patient variability (n = %d)', nPat)); grid on;
subplot(2,2,4);
plot(PlimS, TinjS/60, 'o-', 'LineWidth', 1.5); grid on;
xlabel('allowed IOP limit (mmHg)'); ylabel('injection time (min)');
title('Safety vs procedure-time trade-off (feedback)');
saveFig(fig, figDir, 'F08_iop');

%% ======================================================================
%  S12. SAFETY: EYE-TWITCH DETECTION AND EMERGENCY RETRACTION  (Modules 3, 9)
%  ======================================================================
fprintf('--- S12. Sudden eye movement with the needle inserted (Site A, q3 = %.0f mm) ---\n', qA(3)*mm);
% false-alarm check: 20 s of fixation (no twitch) with tracker noise
rng(P.seed+30);
nS = 20*P.sens.eye_track_fs;
thn = P.sens.eye_track_sd*randn(1,nS+1);            % measured angle, eye still
wn_ = diff(thn)*P.sens.eye_track_fs;                 % finite-difference velocity
ab = abs(wn_) > P.safe.det_thresh | abs(thn(2:end) - thn(1)) > P.safe.det_angle;
fa = sum(ab(2:end) & ab(1:end-1));
fprintf('False alarms during 20 s of steady fixation: %d (velocity noise sigma = %.1f deg/s, threshold %.0f deg/s)\n', ...
    fa, rad2deg(std(wn_)), rad2deg(P.safe.det_thresh));
modes = {'none','retract','spring'}; mname = {'no safety','retraction','spring release'};
RN = cell(1,3);
for m = 1:3, RN{m} = fp_twitch_sim(P.safe.twitch_T, qA, P, modes{m}, P.seed+31); end
fprintf('Nominal twitch: %.0f deg in %.0f ms -> entry point moves %.2f mm\n', rad2deg(P.safe.twitch_amp), ...
    P.safe.twitch_T*1e3, RN{1}.lat(end)*mm);
fprintf('%-24s %-14s %-14s %-22s %-18s\n','Mode','detect (ms)','exit (ms)','lateral at exit (mm)','min lens clr (mm)');
for m = 1:3
    fprintf('%-24s %-14.1f %-14.1f %-22.3f %-18.3f\n', mname{m}, RN{m}.tdet*1e3, RN{m}.texit*1e3, RN{m}.lat_exit*mm, RN{m}.lens_clr_min*mm);
end
Tsw = [0.025 0.05 0.075 0.1 0.15 0.2 0.3 0.4]; latS = zeros(3,numel(Tsw));
for i = 1:numel(Tsw)
    for m = 1:3
        Ri = fp_twitch_sim(Tsw(i), qA, P, modes{m}, P.seed+40+i); latS(m,i) = Ri.lat_exit;
    end
end
fprintf('Twitch-duration sweep (ms): %s\n', num2str(Tsw*1e3));
for m = 1:3, fprintf('  %-22s lateral (mm): %s\n', mname{m}, num2str(latS(m,:)*mm, '%6.3f')); end
okSpring = all(latS(3,Tsw >= 0.1) <= P.safe.max_lateral);
fprintf('Spring quick-release keeps lateral displacement <= %.1f mm for twitches >= 100 ms: %s\n', ...
    P.safe.max_lateral*mm, pf(okSpring));
fprintf('Depth limits on q3: software limit %.0f mm < mechanical hard stop %.0f mm < slide travel %.0f mm; safe workspace clipped by the virtual fixture (S6)\n\n', ...
    P.depth_soft*mm, P.depth_hardstop*mm, P.qmax(3)*mm);
SUM(end+1,:) = {'Eye-twitch false alarms (20 s fixation)', '0', sprintf('%d', fa), fa == 0};
SUM(end+1,:) = {'Scleral lateral displacement (150 ms twitch, spring)', '<= 0.5 mm', sprintf('%.2f mm', RN{3}.lat_exit*mm), RN{3}.lat_exit <= P.safe.max_lateral};

fig = figure('Color','w','Position',[50 50 1150 420]);
subplot(1,2,1); hold on;
lwN = [4 1.5 1.5];                         % 'no safety' is drawn wide: it lies under 'retraction'
for m = 1:3, plot(RN{m}.t*1e3, RN{m}.lat*mm, 'LineWidth', lwN(m)); end
plot(RN{1}.t*1e3, rad2deg(RN{1}.th)/20, 'k:');
yline(P.safe.max_lateral*mm, 'r--', 'limit'); xline(RN{2}.tdet*1e3, 'k--', 'detected');
grid on; xlim([0 400]); xlabel('t (ms)'); ylabel('lateral scleral displacement (mm)');
legend([mname, {'eye angle / 20 (deg)'}], 'Location', 'southeast'); title('Nominal twitch (150 ms, 5 deg)');
subplot(1,2,2);
plot(Tsw*1e3, latS*mm, 'o-', 'LineWidth', 1.5); yline(P.safe.max_lateral*mm, 'r--', 'limit');
grid on; xlabel('twitch duration (ms)'); ylabel('lateral displacement at needle exit (mm)');
legend(mname, 'Location', 'southwest'); title('Safety layer vs speed of the eye movement');
saveFig(fig, figDir, 'F09_safety_twitch');

%% ======================================================================
%  S13. RISK ASSESSMENT, STOPPING DISTANCE, ETHICS  (Module 9)
%  hazard -> risk (S x L) -> control -> verification -> responsibility
%  ======================================================================
fprintf('--- S13. Module 9 safety case: stopping distance, risk register, ethics ---\n');
% (a) stopping distance  S_stop = v*Tr + v^2/(2a)   (M9 slides 10-12)
stp = struct('name', {'Arm transport (fast)', 'Arm docking (slow)', 'Needle insertion (q3)'}, ...
    'v',  {P.arm.vmax, P.arm.vdock, P.vmax(3)}, ...
    'Tr', {P.safe.arm_Tr, P.safe.arm_Tr, P.safe.latency}, ...
    'a',  {P.safe.arm_abrake, P.safe.arm_abrake, P.amax(3)}, ...
    'mrg',{P.safe.arm_margin, P.safe.ins_margin, P.safe.ins_margin}, ...
    'avail', {P.safe.predock, P.safe.predock, P.depth_hardstop - max(qA(3), qB(3))}, ...
    'availTxt', {'pre-dock standoff', 'pre-dock standoff', 'planned depth -> hard stop'});
fprintf('%-24s %-9s %-7s %-9s %-11s %-11s %-11s %-12s %-28s %s\n', 'Motion', 'v (mm/s)', 'Tr (ms)', 'a (m/s2)', ...
    'vTr (mm)', 'v2/2a (mm)', 'S+mrg (mm)', 'avail (mm)', 'available separation', 'result');
okStop = true;
for i = 1:numel(stp)
    s = stp(i); dR = s.v*s.Tr; dB = s.v^2/(2*s.a); stp(i).dR = dR; stp(i).dB = dB;
    stp(i).req = dR + dB + s.mrg; oki = stp(i).req < s.avail; okStop = okStop && oki;
    fprintf('%-24s %-9.1f %-7.0f %-9.3f %-11.3f %-11.3f %-11.3f %-12.1f %-28s %s\n', s.name, s.v*mm, s.Tr*1e3, s.a, ...
        dR*mm, dB*mm, stp(i).req*mm, s.avail*mm, s.availTxt, pf(oki));
end
fprintf('Arm transport speed %.0f mm/s is below the 250 mm/s reduced speed of ISO 10218-1; the insertion drive also has a mechanical hard stop (%.0f mm).\n', ...
    P.arm.vmax*mm, P.depth_hardstop*mm);
fprintf('Simplified estimate only: a certified design needs measured stopping tests (M9 slide 11).\n\n');
SUM(end+1,:) = {'Protective-stop distance < available separation', 'all 3 motions', pf(okStop), okStop};

% (b) protective stop vs emergency stop (M9 slide 9)
fprintf('Stop functions:\n');
fprintf('  Protective stop (automatic): eye-twitch detector (S12), force-sensor puncture/overload, IOP > %d mmHg alarm,\n', P.iop.Palarm);
fprintf('    stereo-registration residual check before docking, joint/workspace limits (virtual fixture, S6).\n');
fprintf('  Emergency stop (manual): red E-stop on the arm cart + surgeon foot pedal; removes drive power,\n');
fprintf('    arm brakes engage, insertion lead screw self-locks (eta_back = %.2f) - it does not replace the protective stops.\n\n', etaBack);

% (c) risk register with severity x likelihood before / after controls (M9 slides 5-6, 21, 23)
p95r = prctile95(eReg);
RR = {
% id hazard                         affected               S Lb La  control                                          verification (this run)
 1, 'Needle advanced too deep',     'patient (retina)',    5, 3, 1, 'mechanical hard stop + software limit + self-locking screw', ...
    sprintf('stop %.2f mm < %.1f mm to hard stop; eta_back %.2f', stp(3).req*mm, stp(3).avail*mm, etaBack), 'hard-stop setting must be inspected before each case'
 2, 'Eye twitch, needle inside',    'patient (sclera, lens)', 4, 3, 1, 'eye tracker 250 Hz + spring quick-release retraction', ...
    sprintf('lateral %.2f mm <= %.1f mm; %d false alarms / 20 s', RN{3}.lat_exit*mm, P.safe.max_lateral*mm, fa), 'twitches < 100 ms exceed the limit (see F09)'
 3, 'Arm hits the face in transport','patient, surgeon, nurse', 4, 3, 1, 'RRT-Connect path (25 mm) + speed & separation monitoring + stop', ...
    sprintf('clearance %.1f mm, stop margin %.1f mm (S7b); S_stop %.1f < %.0f mm', armClr, armSSM, stp(1).req*mm, stp(1).avail*mm), sprintf('braking %.1f m/s^2 from the drive model (S9b); not yet measured', aCapMin)
 4, 'Eye located wrongly (camera)', 'patient (entry site)', 4, 2, 1, 'stereo residual check, refuse docking if residual is high', ...
    sprintf('entry error 95th pct %.3f mm <= %.2f mm', p95r*mm, P.cam.req*mm), 'validated only in simulation (lighting, eyelid occlusion)'
 5, 'IOP spike during injection',   'patient (optic nerve)', 4, 3, 1, 'pressure feedback flow control + 40 mmHg alarm', ...
    sprintf('peak IOP %.1f mmHg < 25 mmHg', RI{3}.Ppeak), 'outflow facility varies between patients'
 6, 'Needle touches lens in re-orientation','patient (lens)', 5, 2, 1, 'C-space planning + dense collision check', ...
    sprintf('min clearance %.2f mm (%s)', clrDense*mm, pf(denseOK)), 'eye model is nominal (no patient-specific scan)'
 7, 'RCM drift stresses the sclera', 'patient (sclera)',  3, 3, 2, 'pivot calibration + tightened tolerances', ...
    sprintf('nominal tol. %.1f um FAIL; %.2fx tol. %.1f um', prctile95(dCal)*um, tolOK, p95tol*um), 'needs the tighter tolerance class -> cost'
 8, 'Power loss with needle inside','patient',             5, 2, 1, 'self-locking screw + arm brakes + manual retract knob', ...
    sprintf('eta_back = %.2f < 0', etaBack), 'manual retraction procedure needs training'
 9, 'Spring / gravity energy in maintenance','maintenance technician', 3, 3, 1, 'lockout: latch pin, release spring, support arm, verify zero energy', ...
    'procedure + inspection checklist', 'administrative control only - not tested here'};
bandNames = {'low','medium','high','critical'}; band = @(r) bandNames{1 + (r > 4) + (r > 9) + (r > 16)};
fprintf('Risk register (score = S x L; bands 1-4 low, 5-9 medium, 10-16 high, 17-25 critical):\n');
fprintf('%-2s %-38s %-24s %-15s %-15s %-62s %s\n', '#', 'Hazard', 'Affected party', 'before (S,L=R)', 'after (S,L=R)', 'Control', 'Verification evidence');
for i = 1:size(RR,1)
    S_ = RR{i,4}; rb = S_*RR{i,5}; ra = S_*RR{i,6};
    fprintf('%-2d %-38s %-24s %d,%d=%-2d %-8s %d,%d=%-2d %-8s %-62s %s\n', RR{i,1}, RR{i,2}, RR{i,3}, S_, RR{i,5}, rb, band(rb), ...
        S_, RR{i,6}, ra, band(ra), RR{i,7}, RR{i,8});
end
fprintf('Residual limitations:\n');
for i = 1:size(RR,1), fprintf('  #%d %s: %s\n', RR{i,1}, RR{i,2}, RR{i,9}); end
okRes = all([RR{:,4}].*[RR{:,6}] <= 9);
SUM(end+1,:) = {'Residual risk of all register hazards', '<= medium (9)', sprintf('max %d', max([RR{:,4}].*[RR{:,6}])), okRes};

% (d) ethics and professional responsibility (M9 slides 16-20)
fprintf('\nEthical / professional analysis (decision - people - evidence - responsibility):\n');
fprintf('  Decision    : allow robot-assisted injection only inside the validated limits above, with the surgeon in the loop\n');
fprintf('                (teleoperation, S6); any failed check -> stop, verify or redesign, never bypass a safeguard.\n');
fprintf('  People      : patient, surgeon, nurse/assistant, maintenance technician, hospital (owner of the device).\n');
fprintf('  Evidence    : S14 summary table + this register (simulation tests, calculations, cited limits).\n');
fprintf('  Responsibility: the surgeon approves each case; the design engineer owns the safety case and residual risks;\n');
fprintf('                the hospital owns maintenance and lockout procedures. Clinical use needs IEC 60601 / ISO 14971 approval.\n');
fprintf('  Perception  : the stereo camera is not a safety proof (like the AI case, slide 19) - docking is refused when the\n');
fprintf('                residual is large, and the surgeon confirms the entry point on the image.\n');
fprintf('  Privacy     : the camera sees the patient face - images are processed on-line, only the eye pose is logged.\n');
fprintf('  Honesty     : the RCM drift FAIL at nominal tolerances is reported, not hidden; all results are simulation only.\n\n');

fig = figure('Color','w','Position',[50 50 1250 520]);
subplot(1,2,1); hold on; axis equal; box on;
[Lg, Sg] = meshgrid(1:5, 1:5); Rg = Lg.*Sg;
cmap = [0.55 0.80 0.55; 1.00 0.90 0.45; 0.98 0.65 0.35; 0.90 0.35 0.35];
for l = 1:5
    for s = 1:5
        r = l*s; c = cmap(1 + (r > 4) + (r > 9) + (r > 16), :);
        rectangle('Position', [l-0.5 s-0.5 1 1], 'FaceColor', c, 'EdgeColor', 'w', 'LineWidth', 1.5);
        text(l + 0.32, s + 0.32, num2str(r), 'FontSize', 7, 'Color', [0.3 0.3 0.3], 'HorizontalAlignment', 'center');
    end
end
Sall = [RR{:,4}];
for i = 1:size(RR,1)
    same = find(Sall == RR{i,4}); rk = find(same == i);          % spread hazards of equal severity vertically
    xb = RR{i,5}; xa = RR{i,6}; y = RR{i,4} + (rk - (numel(same)+1)/2)*0.21;
    quiver(xb, y, xa - xb, 0, 0, 'k', 'LineWidth', 1, 'MaxHeadSize', 0.4);
    plot(xb, y, 'ko', 'MarkerFaceColor', 'w', 'MarkerSize', 11); text(xb, y, num2str(i), 'FontSize', 7, 'HorizontalAlignment', 'center');
    plot(xa, y, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 11); text(xa, y, num2str(i), 'FontSize', 7, 'Color', 'w', 'HorizontalAlignment', 'center');
end
xlim([0.5 5.5]); ylim([0.5 5.5]); xticks(1:5); yticks(1:5);
xlabel('likelihood L (1 rare - 5 frequent)'); ylabel('severity S (1 negligible - 5 catastrophic)');
title({'Risk register: before (white) -> after controls (black)', 'green low | yellow medium | orange high | red critical'});
subplot(1,2,2);
Y = [[stp.dR]; [stp.dB]; [stp.mrg]].'*mm; Av = [stp.avail]*mm;
hb = barh(1:3, Y, 'stacked'); hold on; set(gca, 'XScale', 'log');
hb(1).FaceColor = [0.3 0.5 0.8]; hb(2).FaceColor = [0.9 0.5 0.2]; hb(3).FaceColor = [0.7 0.7 0.7];
plot(Av, 1:3, 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
for i = 1:3, text(Av(i)*1.35, i, sprintf('%.1f mm available', Av(i)), 'Color', 'r', 'FontSize', 8); end
yticks(1:3); yticklabels({stp.name}); grid on; xlim([0.01 600]);
xlabel('distance (mm, log scale)'); legend({'reaction vT_r', 'braking v^2/2a', 'margin', 'separation'}, 'Location', 'southeast');
title('Protective-stop distance  S = vT_r + v^2/(2a) + margin');
saveFig(fig, figDir, 'F12_risk_stopping');

%% ======================================================================
%  S14. SUMMARY AND FMEA  (Module 9)
%  ======================================================================
fprintf('--- S14. Performance summary vs proposal objectives ---\n');
fprintf('%-52s %-20s %-22s %s\n','Metric','Target','Achieved','Result');
for i = 1:size(SUM,1)
    fprintf('%-52s %-20s %-22s %s\n', SUM{i,1}, SUM{i,2}, SUM{i,3}, pf(SUM{i,4}));
end
fprintf('\n--- FMEA (design decisions taken in this project) ---\n');
fmea = {
 'Needle too deep (software fault)', 'Retinal puncture', 'q3 encoder + limit', 'Mechanical hard stop 14 mm + self-locking screw'
 'Lateral force on sclera', 'Scleral tear', 'Kinematic design', 'Mechanical RCM; drift < 10 um after calibration'
 'Hand tremor', 'Tip oscillation', 'Hand-controller sensor', '1:10 scaling + 2nd-order 2 Hz LPF (>90 %)'
 'IOP spike during injection', 'Optic-nerve / CRAO risk', 'Line pressure + model', 'IOP feedback flow control, 40 mmHg alarm'
 'Sudden eye movement', 'Tear / lens touch', 'Eye tracker 250 Hz', 'Detection + spring quick-release retraction'
 'Needle on lens / implant', 'Cataract / implant damage', 'Collision checker', 'C-space planning with clearances'
 'Assembly errors', 'RCM drift', 'Pivot calibration', 'Tolerance allocation from Monte Carlo'
 'Second scleral puncture', 'Endophthalmitis risk', 'Procedure plan', 'Re-orient inside the eye (A -> B)'
 'Eye located wrongly', 'Entry outside pars plana', 'Stereo residual check', 'Camera registration < 0.25 mm + calibration'
 'Arm collides during transport', 'Patient / staff injury', 'Speed + path limits', 'RRT-Connect path (25 mm clearance) + speed and separation monitoring, 120 mm/s cap'};
fprintf('%-34s %-26s %-24s %s\n','Failure mode','Effect','Detection','Mitigation');
for i = 1:size(fmea,1), fprintf('%-34s %-26s %-24s %s\n', fmea{i,:}); end
fprintf('\nTotal runtime: %.1f s\n', toc(tStart));
diary off;
save(fullfile(resDir, 'results.mat'), 'SUM', 'P', 'tipE', 'RI', 'latS', 'Tsw', 'dNom', 'dCal', 'att', 'orders', 'eReg', 'qDockArm');

%% ======================================================================
%  Local functions
%  ======================================================================
function saveFig(fig, figDir, name)
% legends stay inside the axes: small, opaque and placed where they cover the least data
for lg = findobj(fig, 'Type', 'legend').'
    if strcmp(lg.Location, 'none'), continue; end
    if contains(lg.Location, 'outside')
        lg.NumColumns = 1;
        lg.Location = 'best';
    end
    lg.FontSize = 7;
    lg.Color = 'w';
    lg.EdgeColor = [0.55 0.55 0.55];
    lg.LineWidth = 0.5;
    lg.ItemTokenSize = [14 9];
end
drawnow;
print(fig, fullfile(figDir, name), '-dpng', '-r200');
end

function s = pick2(c, a, b)
if c, s = a; else, s = b; end
end

function s = pf(ok)
if ok, s = 'PASS'; else, s = 'FAIL'; end
end

function xd = freePlant(x, P)
[M, c, G] = fp_dynamics(x(1:3), x(4:6), P);
xd = [x(4:6); M \ (-c - G - P.Bv*x(4:6))];
end

function n = marginGradient(p, qprev, P)
% unit gradient (in {0}) of the smallest safety margin at tip position p
h = 1e-6; g = zeros(3,1);
for i = 1:3
    e = zeros(3,1); e(i) = h;
    g(i) = (margin(p+e) - margin(p-e))/(2*h);
end
n = g/max(norm(g), eps);
    function m = margin(pp)
        q = fp_ik(pp, P, qprev);
        [~, cl] = fp_lens_check(q, P);
        m = min([cl(1) - P.eye.safe_lens, cl(2) - P.eye.safe_retina, cl(3) - P.eye.safe_keepout]);
    end
end

function c = pathClearance(pp, d, P)
% minimum lens / keep-out clearance along a C-space polyline (deg)
c = inf;
for k = 1:size(pp,2)-1
    for s = linspace(0,1,10)
        q = pp(:,k) + s*(pp(:,k+1) - pp(:,k));
        [~, cl] = fp_lens_check([deg2rad(q); d], P);
        c = min(c, min(cl([1 3])));
    end
end
end

function [Q,V,A,ph] = addSeg(Q,V,A,ph,kind,q0,qf,opt,dt,id)
if strcmp(kind,'trapezoid'), [~,~,~,T] = fp_traj(kind,q0,qf,0,opt); else, T = opt.T; end
t = 0:dt:ceil(T/dt)*dt;
[q,qd,qdd] = fp_traj(kind,q0,qf,t,opt);
if isempty(Q), Q = q; V = qd; A = qdd; ph = id*ones(1,numel(t));
else, Q = [Q q(:,2:end)]; V = [V qd(:,2:end)]; A = [A qdd(:,2:end)]; ph = [ph id*ones(1,numel(t)-1)]; end
end

function [Q,V,A,ph] = addHold(Q,V,A,ph,q,T,dt,id)
n = round(T/dt);
Q = [Q repmat(q,1,n)]; V = [V zeros(numel(q),n)]; A = [A zeros(numel(q),n)]; ph = [ph id*ones(1,n)];
end

function phaseLines(t, ph)
k = find(diff(ph) ~= 0);
for i = 1:numel(k), xline(t(k(i)), ':', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off'); end
end

function [dNom, dCal] = rcmMonteCarlo(P, sg, nRob, nCfg, seed)
% For every simulated robot: random DH errors, nCfg random poses, needle
% lines (tip, axis). dNom = max distance of the lines to the nominal RCM;
% dCal = max distance to the least-squares pivot point (calibrated RCM).
rng(seed);
dNom = zeros(1,nRob); dCal = dNom;
for r = 1:nRob
    dDH = [sg.th*randn(3,1) sg.d*randn(3,1) sg.a*randn(3,1) sg.al*randn(3,1)];
    A = zeros(3); b = zeros(3,1); L = zeros(6,nCfg);
    for k = 1:nCfg
        q = [deg2rad(-30+60*rand); deg2rad(60+60*rand); 1e-3+9e-3*rand];
        T = fp_fk(q, P, dDH); u = T(1:3,3); p = T(1:3,4);
        L(:,k) = [p; u];
        Pm = eye(3) - u*u.';  A = A + Pm;  b = b + Pm*p;
    end
    xs = A\b;
    dn = zeros(1,nCfg); dc = dn;
    for k = 1:nCfg
        p = L(1:3,k); u = L(4:6,k);
        dn(k) = norm((eye(3)-u*u.')*(0 - p));
        dc(k) = norm((eye(3)-u*u.')*(xs - p));
    end
    dNom(r) = max(dn); dCal(r) = max(dc);
end
end

function v = prctile95(x)
x = sort(x(:)); v = x(max(1, ceil(0.95*numel(x))));
end
