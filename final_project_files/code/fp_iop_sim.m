function R = fp_iop_sim(strategy, P, opt)
%FP_IOP_SIM  Fluidic injection + intraocular-pressure (IOP) dynamics.
%
%   Eye pressure-volume relation (Friedenwald / Pallikaris ocular rigidity):
%       dP/dV = K * P            ->  P = P0 * exp(K * dV)   (closed eye)
%   Aqueous outflow (Goldmann):  Qout = C (P - Pev) + Fu,  production
%   Fprod = C (P0 - Pev) + Fu keeps the eye at P0 before injection.
%   =>  dP/dt = K * P * (Qin + Fprod - Qout)   [P mmHg, Q uL/s, K 1/uL]
%
%   Needle hydraulics (Hagen-Poiseuille, laminar):
%       dP_needle = 128 mu L Q / (pi D^4),  P_line = P_eye + dP_needle
%   A pressure transducer in the syringe line measures P_line (noisy);
%   the controller ESTIMATES the eye pressure as P_line - R_needle*Q
%   (model-based sensor fusion, Module 3).
%
%   strategy: 'bolus'    - whole dose in opt.Tbolus seconds (manual-like)
%             'constant' - whole dose at constant rate over opt.Tconst s
%             'feedback' - IOP-regulated flow: Q = Qout_est + Kp (Pset - P_est),
%                          clipped to [0, opt.Qcap]; Pset = Plim - margin
%   opt.Tend : total simulated time (s);  opt.K, opt.P0 override P.iop
I = P.iop;
if isfield(opt,'K'),  I.K  = opt.K;  end
if isfield(opt,'P0'), I.P0 = opt.P0; end
if ~isfield(opt,'noise'), opt.noise = true; end
mmHg = 133.322;                                   % Pa per mmHg
Rn = 128*I.mu*P.needle.L/(pi*P.needle.ID^4);      % Pa s / m^3
Rn_mmHg_per_uLs = Rn*1e-9/mmHg;                   % mmHg per (uL/s)

dt = 0.01; t = 0:dt:opt.Tend; N = numel(t);
Pe = zeros(1,N); Q = zeros(1,N); V = zeros(1,N); Pl = zeros(1,N); Pest = zeros(1,N);
Pe(1) = I.P0; Vin = 0; done = false; Tinj = NaN;
C = I.C/60; Fu = I.Fu/60;                          % uL/(s mmHg), uL/s
Pset = I.Plim - opt.margin;
for k = 1:N
    % ---- flow command ----
    if done
        q = 0;
    else
        switch strategy
            case 'bolus',    q = I.Vdose/opt.Tbolus;
            case 'constant', q = I.Vdose/opt.Tconst;
            case 'feedback'
                qprev = 0; if k > 1, qprev = Q(k-1); end
                plm = Pe(k) + Rn_mmHg_per_uLs*qprev + opt.noise*P.sens.iop_sd*randn;
                pe_hat = plm - Rn_mmHg_per_uLs*qprev;     % model-based estimate
                if k > 1, pe_hat = 0.9*Pest(k-1) + 0.1*pe_hat; end  % 10 Hz smoothing
                Pest(k) = pe_hat;
                qout_hat = C*max(pe_hat - I.Pev,0) + Fu;
                q = qout_hat + opt.Kp*(Pset - pe_hat);
                q = min(max(q,0), opt.Qcap);
        end
        if Vin + q*dt >= I.Vdose
            q = (I.Vdose - Vin)/dt; done = true; Tinj = t(k) + dt;
        end
    end
    Q(k) = q; Pl(k) = Pe(k) + Rn_mmHg_per_uLs*q;
    if ~strcmp(strategy,'feedback'), Pest(k) = Pe(k); end
    Vin = Vin + q*dt; V(k) = Vin;
    % ---- eye dynamics (RK2) ----
    if k < N
        % aqueous production balances outflow at baseline P0 -> dP/dt = 0
        % at P0 when no drug is injected
        prod = C*(I.P0 - I.Pev) + Fu;
        f = @(p) I.K*p*(q + prod - (C*max(p - I.Pev,0) + Fu));
        k1 = f(Pe(k)); k2 = f(Pe(k) + dt*k1);
        Pe(k+1) = Pe(k) + dt/2*(k1 + k2);
    end
end
R.t = t; R.P = Pe; R.Q = Q; R.V = V; R.Pline = Pl; R.Pest = Pest;
R.Tinj = Tinj; R.Ppeak = max(Pe); R.Rn = Rn; R.Rn_mmHg = Rn_mmHg_per_uLs;
[~, ip] = max(Pe); R.tpeak = t(ip);
% time to return within 2 mmHg of baseline after the peak
kk = find(t > R.tpeak & Pe <= I.P0 + 2, 1);
if isempty(kk), R.Trecover = NaN; else, R.Trecover = t(kk); end
R.time_above_lim = sum(Pe > I.Plim)*dt;
R.dPmax_needle = Rn_mmHg_per_uLs*max(Q);
A_syr = pi*P.iop.syr_D^2/4;
R.Fplunger = max(Pl)*mmHg*A_syr;                  % N, plunger force (gauge)
R.Re = 4*1000*max(Q)*1e-9/(pi*P.needle.ID*I.mu);  % Reynolds number (rho=1000)
end
