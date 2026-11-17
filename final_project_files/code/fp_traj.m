function [q, qd, qdd, T] = fp_traj(kind, q0, qf, t, opt)
%FP_TRAJ  Point-to-point trajectories (Module 7), applied element-wise.
%   kind = 'quintic'   : rest-to-rest 5th-order polynomial of duration opt.T
%   kind = 'trapezoid' : trapezoidal (or triangular) velocity profile using
%                        limits opt.vmax, opt.amax; all joints synchronised
%                        to the slowest joint.
%   t : row vector of time samples (s). Outputs are n x numel(t).
%   T : total duration actually used.
q0 = q0(:); qf = qf(:); t = t(:).';
n = numel(q0); D = qf - q0;
switch kind
    case 'quintic'
        T = opt.T;
        s = min(max(t/T,0),1);
        h   = 10*s.^3 - 15*s.^4 + 6*s.^5;
        hd  = (30*s.^2 - 60*s.^3 + 30*s.^4)/T;
        hdd = (60*s - 180*s.^2 + 120*s.^3)/T^2;
        hd(t > T) = 0; hdd(t > T) = 0;
        q = q0 + D.*h; qd = D.*hd; qdd = D.*hdd;
    case 'trapezoid'
        vmax = opt.vmax(:); amax = opt.amax(:);
        % per-joint minimum time, then synchronise to the slowest joint
        Tj = zeros(n,1);
        for j = 1:n
            d = abs(D(j));
            if d < eps, Tj(j) = 0; continue; end
            if d > vmax(j)^2/amax(j)
                Tj(j) = d/vmax(j) + vmax(j)/amax(j);   % trapezoid
            else
                Tj(j) = 2*sqrt(d/amax(j));              % triangle
            end
        end
        T = max(Tj);
        if isfield(opt,'T') && opt.T > T, T = opt.T; end
        q = zeros(n,numel(t)); qd = q; qdd = q;
        for j = 1:n
            d = abs(D(j)); sg = sign(D(j));
            if d < eps, q(j,:) = q0(j); continue; end
            % choose cruise speed so that the move takes exactly T with amax(j)
            disc = (amax(j)*T)^2 - 4*amax(j)*d;
            v = (amax(j)*T - sqrt(max(disc,0)))/2;
            ta = v/amax(j); a = amax(j);
            for k = 1:numel(t)
                tk = min(max(t(k),0),T);
                if tk < ta
                    s = 0.5*a*tk^2; sd = a*tk; sdd = a;
                elseif tk < T - ta
                    s = 0.5*a*ta^2 + v*(tk-ta); sd = v; sdd = 0;
                else
                    tr = T - tk;
                    s = d - 0.5*a*tr^2; sd = a*tr; sdd = -a;
                end
                if t(k) >= T, sd = 0; sdd = 0; s = d; end
                q(j,k) = q0(j) + sg*s; qd(j,k) = sg*sd; qdd(j,k) = sg*sdd;
            end
        end
    otherwise
        error('unknown trajectory kind');
end
end
