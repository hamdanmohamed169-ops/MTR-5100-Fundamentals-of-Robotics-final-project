function S = fp_simulate(t, Qd, Vd, Ad, P, ctrl)
%FP_SIMULATE  Closed-loop joint-space simulation of the RCM robot (Module 7).
%   t          : 1xN time vector, uniform step
%   Qd,Vd,Ad   : 3xN desired joint position / velocity / acceleration
%   ctrl.type  : 'PID' (decentralised, model-free) or 'CTC' (computed torque)
%   ctrl.wn, ctrl.zeta : desired closed-loop natural frequency / damping
%   ctrl.model_scale   : mass scale used inside the CTC model (1 = exact)
%   ctrl.needle        : true -> needle/tissue force acts on joint 3 (insertion)
%
%   Plant  : M(q)qdd + c + G + Bv qd = sat(tau) - J'F_needle   (RK4, ZOH tau)
%   Sensor : motor encoders quantise q (resolution 2*pi/(cpr*Nr));
%            velocity = filtered backward difference of the quantised q.
%   Actuator saturation: peak motor torque = 3 x rated, through N and eta.
N  = numel(t); dt = t(2) - t(1); n = size(P.DH,1);
rv = P.DH(:,5) == 0;                    % revolute joints (angle wrapping)
if ~isfield(ctrl,'model_scale'), ctrl.model_scale = 1; end
if ~isfield(ctrl,'needle'), ctrl.needle = true; end

% --- actuator limits at the joint side (Module 2) ---
A = P.act;
tau_max = [3*A.Tm_rated.*A.N.*A.eta, ...
           3*A.Tm4_rated*2*pi*A.eta_s*A.N4/A.lead].';

% --- encoder resolution at the joint side (Module 3) ---
res = 2*pi./(P.sens.enc_cpr*A.Nr(:));

% --- gains ---
wn = ctrl.wn; ze = ctrl.zeta;
M0 = fp_dynamics(Qd(:,1), zeros(n,1), P, ctrl.model_scale);
Kp = wn^2*ones(n,1); Kd = 2*ze*wn*ones(n,1); Ki = zeros(n,1);
if strcmp(ctrl.type,'PID')
    m = diag(M0);
    Kp = m*wn^2;  Kd = m*2*ze*wn;  Ki = Kp*wn/5;
end

% --- logs ---
S.t = t; S.q = zeros(n,N); S.qd = zeros(n,N); S.tau = zeros(n,N);
S.tau_cmd = zeros(n,N); S.F = zeros(1,N); S.punct = false(1,N);
S.tau_max = tau_max; S.res = res;

x  = [Qd(:,1); Vd(:,1)];
qm_prev = round(x(1:n)./res).*res; vf = zeros(n,1);
ei = zeros(n,1);
nst.punctured = false;
alpha_v = exp(-2*pi*50*dt);            % 50 Hz velocity filter
for k = 1:N
    q = x(1:n); qd = x(n+1:2*n);
    % ---------- measurement ----------
    qm = round(q./res).*res;
    vraw = (qm - qm_prev)/dt; qm_prev = qm;
    vf = alpha_v*vf + (1-alpha_v)*vraw;
    if k == 1, vf = Vd(:,1); end
    e  = Qd(:,k) - qm;  ed = Vd(:,k) - vf;
    e(rv) = atan2(sin(e(rv)), cos(e(rv)));
    % ---------- controller ----------
    switch ctrl.type
        case 'PID'
            ei = ei + e*dt;
            tc = Kp.*e + Kd.*ed + Ki.*ei;
        case 'CTC'
            [Mh, ch, Gh] = fp_dynamics(qm, vf, P, ctrl.model_scale);
            tc = Mh*(Ad(:,k) + Kd.*ed + Kp.*e) + ch + Gh + P.Bv*vf;
    end
    tau = min(max(tc, -tau_max), tau_max);
    % ---------- needle force (held over the step) ----------
    F = 0;
    if ctrl.needle
        [F, nst] = fp_needle_force(q(n), qd(n), nst, P);
    end
    tau_ext = [zeros(n-1,1); -F];
    % ---------- log ----------
    S.q(:,k) = q; S.qd(:,k) = qd; S.tau(:,k) = tau; S.tau_cmd(:,k) = tc;
    S.F(k) = F; S.punct(k) = nst.punctured;
    % ---------- integrate (RK4) ----------
    if k < N
        f = @(xx) plant(xx, tau + tau_ext, P);
        k1 = f(x); k2 = f(x + dt/2*k1); k3 = f(x + dt/2*k2); k4 = f(x + dt*k3);
        x = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    end
end
S.e = Qd - S.q;
S.e(rv,:) = atan2(sin(S.e(rv,:)), cos(S.e(rv,:)));
end

function xd = plant(x, tau, P)
n = size(P.DH,1); q = x(1:n); qd = x(n+1:2*n);
[M, c, G] = fp_dynamics(q, qd, P);
qdd = M \ (tau - c - G - P.Bv*qd);
xd = [qd; qdd];
end
