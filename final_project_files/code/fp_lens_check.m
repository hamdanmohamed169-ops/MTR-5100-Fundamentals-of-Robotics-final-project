function [ok, clr, info] = fp_lens_check(q, P, thetaEye)
%FP_LENS_CHECK  Collision / clearance test of the needle inside the eye
%   (Module 8 collision checking). Works in the eye frame {E}.
%   q        : joint vector [q1 yaw; q2 pitch; q3 depth]
%   thetaEye : optional eye rotation about xE (rad) for the twitch study
%   ok   : true if the inserted shaft clears the lens by safe_lens AND the
%          tip keeps safe_retina from the retina
%          AND the shaft keeps safe_keepout from the keep-out sphere
%   clr  : [lens, retina, keep-out] clearances in metres (lens clearance uses
%          a radial approximation of point-to-ellipsoid distance)
if nargin < 3, thetaEye = 0; end
[p0, ~] = fp_fk_geom(q);
R = P.T_E0(1:3,1:3); o = P.T_E0(1:3,4);
a = o;                    % scleral entry (RCM) in {E}
b = o + R*p0;             % needle tip in {E}
if thetaEye ~= 0          % express needle in the rotated eye's own frame
    c = cos(thetaEye); s = sin(thetaEye);
    Rx = [1 0 0; 0 c -s; 0 s c];
    a = Rx.'*a; b = Rx.'*b;
end
if q(3) <= 0              % tip still outside the globe
    ok = true; clr = [inf inf inf]; info.tipE = b; return
end
ax = P.eye.lens_ax; cL = P.eye.lens_c;
ns = 25; s = linspace(0,1,ns);
pts = a + (b - a).*s;
dl = inf; minScaled = inf;
for k = 1:ns
    v = pts(:,k) - cL;
    sc = norm(v./ax);
    minScaled = min(minScaled, sc);
    dl = min(dl, norm(v)*(1 - 1/sc));    % radial distance to the surface
end
% exact segment test against the inflated ellipsoid (unit sphere in scaled space)
axi = ax + P.eye.safe_lens;
A = (a - cL)./axi; B = (b - cL)./axi; d = B - A;
tt = min(max(-(A.'*d)/(d.'*d), 0), 1);
hitInflated = norm(A + tt*d) < 1;
dr = P.eye.Rret - norm(b);           % retina clearance of the tip
% keep-out sphere: exact segment-to-point distance
cK = P.eye.keepout_c; ab = b - a;
tk = min(max(((cK - a).'*ab)/(ab.'*ab), 0), 1);
dk = norm(a + tk*ab - cK) - P.eye.keepout_r;
ok = ~hitInflated && dr >= P.eye.safe_retina && dk >= P.eye.safe_keepout;
clr = [dl dr dk];
info.tipE = b; info.minScaled = minScaled;
end
