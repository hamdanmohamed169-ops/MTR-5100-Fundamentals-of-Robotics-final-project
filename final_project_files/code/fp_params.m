function P = fp_params()
%FP_PARAMS  All parameters and assumptions of the final project in one place.
%   3-DOF RCM ophthalmic injection robot (R-R-P) carried by a 4-DOF
%   positioning arm, with a stereo camera that locates the eye. SI units unless the
%   field name says otherwise (_deg, _mm, _mmHg, _uL are clinical display
%   units and are converted explicitly where they are used).
%
%   MTR 5100 - Fundamentals of Robotics, Final Project
%   Mohamed Abdelaziz Hamdan

%% ---------------- General ----------------------------------------------
P.g     = 9.81;                 % m/s^2
P.seed  = 5100;                 % RNG seed -> reproducible results
P.dt    = 1e-3;                 % s, control/simulation step (1 kHz)

%% ---------------- Eye model (frame {E}) --------------------------------
% {E}: origin at the eye-globe centre, zE along the visual axis toward the
% cornea, xE temporal, patient supine and looking up -> world "up" = +zE.
P.eye.R        = 12.0e-3;       % m, globe outer radius (axial length 24 mm)
P.eye.Rret     = 11.0e-3;       % m, inner retina surface radius
P.eye.limbus_r = 5.75e-3;       % m, limbus radius (cornea diameter 11.5 mm)
P.eye.pp_arc   = 4.0e-3;        % m, pars plana entry 4 mm behind limbus (along surface)
P.eye.lens_c   = [0; 0; 6.5e-3];        % m, crystalline lens centre in {E}
P.eye.lens_ax  = [4.6e-3; 4.6e-3; 2.0e-3]; % m, lens semi-axes (9.2 mm dia, 4.0 mm thick)
P.eye.safe_lens   = 1.0e-3;     % m, required clearance needle <-> lens
P.eye.safe_retina = 2.0e-3;     % m, required clearance tip <-> retina
P.eye.gE       = [0; 0; -P.g];  % gravity expressed in {E} (supine patient)

% Pars plana entry point (temporal side, phi = 0)
thL  = asin(P.eye.limbus_r/P.eye.R);          % polar angle of limbus
thE  = thL + P.eye.pp_arc/P.eye.R;            % polar angle of entry point
P.eye.theta_entry = thE;
P.eye.pE = P.eye.R*[sin(thE); 0; cos(thE)];   % entry point (= RCM) in {E}

%% ---------------- Robot base frame {0} at the RCM ----------------------
% x0 points from the entry point to the eye centre (nominal insertion axis),
% z0 is the component of world-up perpendicular to x0, y0 = z0 x x0.
x0 = -P.eye.pE/norm(P.eye.pE);
up = [0;0;1];
z0 = up - (up.'*x0)*x0;  z0 = z0/norm(z0);
y0 = cross(z0, x0);
P.T_E0 = [x0 y0 z0 P.eye.pE; 0 0 0 1];        % pose of {0} in {E}

% Keep-out zone from pre-operative imaging (e.g. a free-floating
% intravitreal implant or vitreous opacity): sphere the needle shaft must
% not touch while it is re-oriented inside the eye (Module 8 obstacle).
pK0 = fp_fk_geom([0; deg2rad(100); 4.5e-3]);
P.eye.keepout_c = P.eye.pE + P.T_E0(1:3,1:3)*pK0;   % m, centre in {E}
P.eye.keepout_r = 0.5e-3;                            % m, radius
P.eye.safe_keepout = 0.2e-3;                         % m, clearance

%% ---------------- Kinematics (standard DH) -----------------------------
% Row i: [theta_offset  d  a  alpha  type]  (type 0 = revolute, 1 = prismatic)
% Spherical R-R-P robot: both revolute axes intersect at the RCM -> needle
% direction u = [c1 s2, s1 s2, c2] in {0}; tip position p = q3*u.
P.DH = [ 0   0   0  -pi/2  0 ;     % joint 1 yaw   (about z0)
         0   0   0   pi/2  0 ;     % joint 2 pitch
         0   0   0   0     1 ];    % joint 3 insertion (d3 = q3)
P.qmin = [deg2rad(-40); deg2rad(50);  -15e-3];
P.qmax = [deg2rad( 40); deg2rad(130);  20e-3];
P.vmax = [deg2rad(20); deg2rad(20); 5e-3];      % rad/s, rad/s, m/s
P.amax = [deg2rad(80); deg2rad(80); 20e-3];     % rad/s^2, rad/s^2, m/s^2
P.q_home = [0; pi/2; -10e-3];        % needle aligned, tip 10 mm outside eye
% Depth limits on q3, ordered so that software always acts before hardware:
%   software limit (13 mm) < mechanical hard stop (14 mm) < slide travel qmax(3) (20 mm)
P.depth_soft     = 13e-3;            % m, software depth limit used by the controller / virtual fixture
P.depth_hardstop = 14e-3;            % m, mechanical hard stop on q3 (inside the slide travel qmax(3))

%% ---------------- Mass properties (for dynamics) -----------------------
% Physical RCM mechanism represented by 3 rigid bodies. COM positions are in
% each body's own DH frame (link i rides on frame i). Masses are typical of
% bench-top ophthalmic robots (total moving mass < 1.5 kg).
P.link(1).m = 0.60;  P.link(1).rc = [0; -0.06; 0];    P.link(1).I = diag([1.2 0.4 1.2])*1e-3;
% pitch arc + tool unit + insertion actuator (one rigid body)
P.link(2).m = 0.60;  P.link(2).rc = [0;  0.0; -0.065]; P.link(2).I = diag([1.05 1.05 0.32])*1e-3;
% Insertion carriage + syringe: its COM is 60 mm behind the needle tip
P.link(3).m = 0.12;  P.link(3).rc = [0;  0;  -0.06];  P.link(3).I = diag([0.05 0.05 0.01])*1e-3;
P.Bv = diag([2e-3 2e-3 2.0]);        % viscous friction (N m s/rad, N m s/rad, N s/m)

%% ---------------- Actuation (Module 2) ---------------------------------
% Joints 1-2: brushless DC + harmonic drive; joint 3: DC motor + lead screw
P.act.N     = [100 100];             % gear ratios joints 1..2
P.act.eta   = [0.75 0.75];           % gearbox efficiency
P.act.Tm_rated = [0.020 0.020];      % N m, rated motor torque (EC-max 16 class)
P.act.wm_rated = [5000 5000]*2*pi/60;% rad/s
P.act.lead  = 0.5e-3;                % m/rev, lead screw (joint 3, insertion)
P.act.eta_s = 0.35;                  % lead-screw efficiency (< 0.5 -> self-locking)
P.act.N4    = 1;                     % motor:screw reduction (direct drive), insertion drive
P.act.Tm4_rated = 0.008;             % N m
P.act.wm4_rated = 6000*2*pi/60;      % rad/s
P.act.SF    = 1.5;                   % design safety factor
P.act.Jm    = [6e-8 6e-8 6e-8];      % kg m^2, rotor inertias (joints 1..3)
% transmission ratio joint->motor: rad/rad for 1-2, rad/m for joint 3
P.act.Nr    = [P.act.N, 2*pi*P.act.N4/P.act.lead];

%% ---------------- Sensing (Module 3) -----------------------------------
P.sens.enc_cpr   = 4096;             % decoded counts/rev at the motor
P.sens.force_sd  = 2e-3;             % N, needle force sensor noise (1 sigma)
P.sens.force_fs  = 1000;             % Hz
P.sens.iop_sd    = 0.3;              % mmHg, pressure transducer noise
P.sens.eye_track_fs = 250;           % Hz, eye tracking (same stereo camera, region-of-interest mode)
P.sens.eye_track_sd = deg2rad(0.01); % rad, eye-angle measurement noise (precision)

%% ---------------- Needle / tissue interaction --------------------------
% 30G needle; piecewise insertion-force model (sclera puncture).
P.needle.OD   = 0.311e-3;   % m
P.needle.ID   = 0.159e-3;   % m
P.needle.L    = 12.7e-3;    % m (1/2 inch)
P.needle.k_scl = 150;       % N/m, scleral wall stiffness before puncture (tenting)
P.needle.F_punct = 0.20;    % N, puncture force (30G, bevelled)
P.needle.f_fric  = 5.0;     % N/m, shaft friction per metre of inserted length
P.needle.c_vit   = 0.5;     % N s/m, vitreous viscous drag

%% ---------------- Tremor and teleoperation -----------------------------
P.trem.amp   = 100e-6;      % m, physiological tremor amplitude at the hand
P.trem.fband = [5 12];      % Hz
P.trem.scale = 0.1;         % motion scaling 1:10
P.trem.fc    = 2.0;         % Hz, low-pass cut-off for the tremor filter

%% ---------------- Fluidics and IOP (Friedenwald/Pallikaris) ------------
P.iop.P0     = 15;          % mmHg, baseline IOP
P.iop.K      = 0.0126;      % 1/uL, ocular rigidity (dP/dV = K*P)
P.iop.C      = 0.30;        % uL/(min mmHg), trabecular outflow facility
P.iop.Pev    = 9;           % mmHg, episcleral venous pressure
P.iop.Fu     = 0.4;         % uL/min, uveoscleral (pressure-independent) outflow
P.iop.Vdose  = 50;          % uL, dose (anti-VEGF standard 0.05 mL)
P.iop.Plim   = 25;          % mmHg, design target from the proposal
P.iop.Palarm = 40;          % mmHg, hard alarm (risk to retinal perfusion)
P.iop.mu     = 5e-3;        % Pa s, drug viscosity (protein solution, ~5 cP)
P.iop.syr_D  = 4.78e-3;     % m, 1 mL syringe barrel inner diameter
P.iop.Qmax   = 100;         % uL/s, pump flow limit

%% ---------------- Safety (Module 9) ------------------------------------
P.safe.twitch_amp  = deg2rad(5);    % rad, sudden eye rotation (micro-saccade/twitch)
P.safe.twitch_T    = 0.150;         % s, nominal twitch duration (sweep 25-400 ms)
P.safe.det_thresh  = deg2rad(30);   % rad/s, eye-velocity threshold for detection
P.safe.det_angle   = deg2rad(1.0);  % rad, eye-displacement threshold (slow drifts)
P.safe.v_retract   = 50e-3;         % m/s, emergency retraction speed (joint 3)
P.safe.a_retract   = 5.0;           % m/s^2, emergency retraction acceleration
P.safe.latency     = 4e-3;          % s, processing + actuation latency after detection
% Design alternative: spring-loaded quick-release retraction (latch + preloaded spring)
P.safe.v_spring    = 0.30;          % m/s, spring retraction speed
P.safe.a_spring    = 60;            % m/s^2, spring retraction acceleration
P.safe.latency_spring = 3e-3;       % s, solenoid latch release time
P.safe.max_lateral = 0.5e-3;        % m, max allowed scleral lateral displacement
% Protective-stop data for the stopping-distance check, M9: S = v*Tr + v^2/(2a) + margin
P.safe.arm_Tr      = 0.050;         % s, arm protective-stop response (monitoring + safety relay + drive)
P.safe.arm_abrake  = 1.0;           % m/s^2, arm controlled-stop deceleration (Cartesian, conservative)
P.safe.arm_margin  = 0.010;         % m, extra safety margin
P.safe.predock     = 0.040;         % m, fast transport ends 40 mm before the entry (demo_live xPre)
P.safe.ins_margin  = 0.5e-3;        % m, extra margin on the insertion axis

%% ---------------- 4-DOF positioning arm (macro stage) -----------------
% Carries the RCM unit from its parking pose to the eye and back. Base yaw
% on a vertical column + a 3-link planar arm (shoulder, elbow, wrist pitch).
% Standard DH, metres:
%   i | theta | d    | a    | alpha
%   1 | q1    | d1   | 0    | +pi/2     base yaw (vertical column)
%   2 | q2    | 0    | a2   | 0         shoulder pitch
%   3 | q3    | 0    | a3   | 0         elbow pitch
%   4 | q4    | 0    | a4   | 0         wrist pitch -> flange frame {4}
% Task space (4 DOF): flange position (x, y, z) + pitch angle phi = q2+q3+q4.
% Closed-form IK: q1 from atan2, planar 2-link IK for (q2, q3) at the
% wrist point, q4 = phi - q2 - q3 -> 4 solutions (front/back x elbow up/down).
P.arm.d1 = 0.400; P.arm.a2 = 0.300; P.arm.a3 = 0.260; P.arm.a4 = 0.120;
P.arm.DH = [0 P.arm.d1 0 pi/2; 0 0 P.arm.a2 0; 0 0 P.arm.a3 0; 0 0 P.arm.a4 0];
P.arm.T_Eb = [-1 0 0 0.250; 0 -1 0 0; 0 0 1 -0.150; 0 0 0 1];   % arm base (table-rail clamp) in {E}, x_b towards the patient
P.arm.dF   = 0.182;          % flange -> RCM point, along the flange x-axis
% RCM unit on the flange: the flange x-axis points down the unit (-z0) and
% the pitch axes are parallel to y0, so the docking pose lies in the arm plane.
x0 = P.T_E0(1:3,1); y0 = P.T_E0(1:3,2); z0 = P.T_E0(1:3,3);
T_E4dock = [-z0 x0 -y0 P.eye.pE + P.arm.dF*z0; 0 0 0 1];
P.arm.T_40 = T_E4dock\P.T_E0;  % fixed transform flange {4} -> RCM base {0}
P.arm.qmin = deg2rad([-170; -40; -170; -180]);
P.arm.qmax = deg2rad([ 170; 180;  170;  180]);
P.arm.vmax = 0.120;          % m/s, Cartesian transport speed limit (well below the 250 mm/s reduced-speed limit of ISO 10218)
P.arm.vdock = 0.010;         % m/s, final docking speed along the insertion axis
P.arm.wmax  = 0.60;          % rad/s, joint speed limit of the arm during transport

% Mass properties of the arm (Modules 2, 6). DH frames {i}: x_i along link i
% (towards the next joint), y_1 up at q = 0, z_i = next joint axis.
% COM rc (m) and inertia I about the COM (kg m^2) are expressed in frame {i}.
%   link 1: turret + shoulder housing (J2 motor + harmonic drive), COM 30 mm below the shoulder axis
%   link 2: upper arm tube 1.0 kg (COM at mid-length) + elbow housing 1.0 kg (J3 drive) at the elbow
%   link 3: forearm tube 0.7 kg + wrist housing 0.7 kg (J4 drive) at the wrist
%   link 4: wrist link + flange
P.arm.link(1) = struct('m', 2.5, 'rc', [0; -0.030; 0],  'I', diag([3.7e-3 2.5e-3 3.7e-3]));
P.arm.link(2) = struct('m', 2.0, 'rc', [-0.075; 0; 0],  'I', diag([1.2e-3 1.88e-2 1.88e-2]));
P.arm.link(3) = struct('m', 1.4, 'rc', [-0.065; 0; 0],  'I', diag([6.0e-4 9.9e-3 9.9e-3]));
P.arm.link(4) = struct('m', 0.4, 'rc', [-0.060; 0; 0],  'I', diag([1.0e-4 5.0e-4 5.0e-4]));
% Carried RCM injection unit (robot links 1.32 kg + J1 housing + quick-connect + cable):
% COM and inertia from the drawn unit geometry in the flange frame {4}
P.arm.payload = struct('m', 1.85, 'rc', [0.1625; -0.0549; -0.0116], ...
    'I', [1.9e-3 1.5e-3 0.2e-3; 1.5e-3 5.6e-3 -0.2e-3; 0.2e-3 -0.2e-3 6.7e-3]);
% Arm drives: brushless motor + harmonic drive per joint = the result of the
% sizing workflow in S9b (smallest catalogue motor / ratio meeting SF 1.5):
%   J1 30 W 1:100, J2 150 W 1:160, J3 70 W 1:160, J4 30 W 1:100
P.arm.act.N        = [100 160 160 100];            % harmonic-drive ratios
P.arm.act.eta      = [0.75 0.75 0.75 0.75];        % gear efficiency
P.arm.act.Tm_rated = [0.05 0.25 0.10 0.05];        % N m, rated (continuous) motor torque
P.arm.act.Tm_peak  = 3*[0.05 0.25 0.10 0.05];      % N m, peak motor torque (3 x rated)
P.arm.act.wm_rated = [6000 4500 5000 6000]*2*pi/60;% rad/s, rated motor speed
P.arm.act.Jm       = [4.0e-6 3.3e-5 1.2e-5 4.0e-6];% kg m^2, rotor inertia
P.arm.act.Tbrake   = [0 58 28 11];                 % N m, power-off brake holding torque at the joint (2 x worst gravity)

%% ---------------- Stereo camera for eye registration (Modules 3-4) -----
% Two calibrated pinhole cameras on a boom above the patient. They detect
% points on the limbus (cornea/sclera border), triangulate them, fit a 3-D
% circle -> eye centre, optical axis and pars plana entry point.
P.cam.C      = [-0.060; 0.110; 0.400];   % m, mid-point of the stereo pair in {E}
P.cam.target = [0; 0; 0.010];             % m, aim point
P.cam.base   = 0.120;                     % m, stereo baseline
P.cam.f      = 3500;                      % px, focal length
P.cam.res    = [1280 1024];               % px, image size
P.cam.sd_px  = 0.25;                      % px, corner/edge detection noise (1 sigma)
P.cam.nPts   = 48;                        % limbus points per image
P.cam.sd_calib_t = 0.05e-3;               % m, camera-to-robot (hand-eye) calibration error
P.cam.sd_calib_r = deg2rad(0.01);         % rad
P.cam.req    = 0.25e-3;                   % m, required entry-point accuracy (pars plana band +-0.5 mm)
end
