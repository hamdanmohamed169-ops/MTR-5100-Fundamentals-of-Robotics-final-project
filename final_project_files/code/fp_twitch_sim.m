function R = fp_twitch_sim(Ttw, q0, P, mode, seed)
%FP_TWITCH_SIM  Sudden eye rotation while the needle is inserted (safety).
%   The eye rotates about xE (worst case: entry point moves perpendicular to
%   the needle) by P.safe.twitch_amp with a minimum-jerk profile of duration
%   Ttw. The robot holds the needle fixed in space unless the safety layer
%   detects the motion:
%     - eye tracker samples at P.sens.eye_track_fs with Gaussian noise,
%     - velocity = difference of consecutive samples,
%     - detection = 2 consecutive samples with velocity > det_thresh OR
%       displacement from the fixation reference > det_angle (slow drifts),
%     - after P.safe.latency, joint 3 (insertion) retracts with a_retract / v_retract.
%   mode = 'none' (no safety layer), 'retract' (motor + lead screw) or
%          'spring' (spring-loaded quick-release, design alternative).
%   Outputs: lateral scleral displacement at the time the needle leaves the
%   eye (or at the end), detection time, exit time, minimum lens clearance.
rng(seed);
dt = 1e-4; t = 0:dt:0.8; N = numel(t);
A = P.safe.twitch_amp;
s = min(t/Ttw,1);
th = A*(10*s.^3 - 15*s.^4 + 6*s.^5);
% eye tracker
Ts = 1/P.sens.eye_track_fs;
ts = 0:Ts:t(end);
thm = interp1(t, th, ts) + P.sens.eye_track_sd*randn(size(ts));
wm = [0 diff(thm)/Ts];
above = abs(wm) > P.safe.det_thresh | abs(thm - thm(1)) > P.safe.det_angle;
kdet = find(above(2:end) & above(1:end-1), 1) + 1;
if isempty(kdet), tdet = inf; else, tdet = ts(kdet); end
% needle line (fixed in space) and entry point
R0 = P.T_E0(1:3,1:3); pE = P.T_E0(1:3,4);
[~, u0] = fp_fk_geom(q0); uE = R0*u0;
q4 = q0(3); v = 0; texit = NaN;
lat = zeros(1,N); q4l = zeros(1,N); clr = inf(1,N);
for k = 1:N
    c = cos(th(k)); sn = sin(th(k));
    Rx = [1 0 0; 0 c -sn; 0 sn c];
    dp = Rx*pE - pE;
    if q4 > 0                         % needle still crossing the sclera
        lat(k) = norm(dp - (dp.'*uE)*uE);
        [~, cl] = fp_lens_check([q0(1:2); q4], P, th(k));
        clr(k) = cl(1);
    else
        lat(k) = lat(max(k-1,1));
        if isnan(texit), texit = t(k); end
    end
    q4l(k) = q4;
    switch mode
        case 'retract'      % motor-driven retraction through the lead screw
            if t(k) >= tdet + P.safe.latency && q4 > -1e-3
                v = max(v - P.safe.a_retract*dt, -P.safe.v_retract);
                q4 = q4 + v*dt;
            end
        case 'spring'       % spring-loaded quick-release retraction
            if t(k) >= tdet + P.safe.latency_spring && q4 > -1e-3
                v = max(v - P.safe.a_spring*dt, -P.safe.v_spring);
                q4 = q4 + v*dt;
            end
    end
end
R.t = t; R.th = th; R.lat = lat; R.q4 = q4l; R.tdet = tdet; R.texit = texit;
R.lat_exit = lat(end); R.lens_clr_min = min(clr);
R.ts = ts; R.wm = wm;
end
