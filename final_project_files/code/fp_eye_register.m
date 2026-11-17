function R = fp_eye_register(P, seed)
%FP_EYE_REGISTER  Locate the eye with the calibrated stereo camera (Modules 3-4).
%   1. Two pinhole cameras (intrinsics K, extrinsics [R|t]) image the limbus
%      (the circular cornea/sclera border); every detected edge point has
%      Gaussian pixel noise P.cam.sd_px.
%   2. Each left/right point pair is triangulated with the linear DLT
%      method (SVD of the 4x4 system A X = 0).
%   3. The camera->robot (hand-eye) calibration error is applied.
%   4. A plane is fitted to the 3-D points (SVD) -> limbus normal n, and a
%      circle is fitted in that plane (algebraic least squares) -> limbus
%      centre c and radius r.
%   5. Eye-model geometry gives the globe centre  o = c - n sqrt(R^2 - rL^2)
%      and the pars plana entry point  pE = o + R (sin(thE) x_t + cos(thE) n),
%      x_t = temporal direction (table lateral axis projected on the plane).
%   6. The estimated robot base frame T_E0_est is built exactly like P.T_E0.
%   R.err_entry, R.err_centre (m) and R.err_axis (rad) compare with the truth.
if nargin < 2, seed = P.seed; end
rng(seed);
C = P.cam;
% ---------------- camera models ----------------
zc = (C.target - C.C)/norm(C.target - C.C);
bx = cross(zc, [0;0;1]); bx = bx/norm(bx);       % baseline direction
K = [C.f 0 C.res(1)/2; 0 C.f C.res(2)/2; 0 0 1];
cams(1).Cc = C.C - C.base/2*bx;                  % left camera centre
cams(2).Cc = C.C + C.base/2*bx;                  % right camera centre
for i = 1:2                                      % both cameras verge on the eye
    zi = (C.target - cams(i).Cc)/norm(C.target - cams(i).Cc);
    xi = cross(zi, [0;0;1]); xi = xi/norm(xi);   % image x (horizontal)
    yi = cross(zi, xi);                          % image y (down in the image)
    Rwc = [xi yi zi];                            % camera axes in {E}
    cams(i).R = Rwc; cams(i).K = K;
    cams(i).Pm = K*[Rwc.' -Rwc.'*cams(i).Cc];    % 3x4 projection matrix
end
% ---------------- true limbus and detections ----------------
rL = P.eye.limbus_r; zL = sqrt(P.eye.R^2 - rL^2);
ang = 2*pi*((0:C.nPts-1) + 0.3*rand(1,C.nPts))/C.nPts;
Xt = [rL*cos(ang); rL*sin(ang); zL*ones(1,C.nPts)];
uv = cell(1,2);
for i = 1:2
    uv{i} = proj(cams(i).Pm, Xt) + C.sd_px*randn(2, C.nPts);
end
% ---------------- triangulation (linear DLT) ----------------
Xs = zeros(3, C.nPts);
for k = 1:C.nPts
    A = zeros(4,4);
    for i = 1:2
        Pm = cams(i).Pm;
        A(2*i-1,:) = uv{i}(1,k)*Pm(3,:) - Pm(1,:);
        A(2*i,:)   = uv{i}(2,k)*Pm(3,:) - Pm(2,:);
    end
    [~,~,V] = svd(A); Xh = V(:,end);
    Xs(:,k) = Xh(1:3)/Xh(4);
end
% ---------------- hand-eye calibration error ----------------
w = randn(3,1); w = w/norm(w);
Rerr = axang(w, C.sd_calib_r*randn); terr = C.sd_calib_t*randn(3,1);
Xm = Rerr*(Xs - C.C) + C.C + terr;              % points as the robot sees them
% ---------------- plane + circle fit ----------------
m = mean(Xm, 2);
[U,~,~] = svd(Xm - m, 'econ');
n = U(:,3); if n.'*(C.C - m) < 0, n = -n; end   % normal towards the camera
e1 = U(:,1); e2 = cross(n, e1);
xy = [e1 e2].'*(Xm - m);
sol = [2*xy.' ones(C.nPts,1)] \ sum(xy.^2, 1).';
cEst = m + e1*sol(1) + e2*sol(2);
rEst = sqrt(sol(3) + sol(1)^2 + sol(2)^2);
% ---------------- eye centre, entry point, robot base frame ----------------
oEst = cEst - n*sqrt(P.eye.R^2 - rL^2);
xt = [1;0;0]; xt = xt - (xt.'*n)*n; xt = xt/norm(xt);
thE = P.eye.theta_entry;
pEst = oEst + P.eye.R*(sin(thE)*xt + cos(thE)*n);
x0 = (oEst - pEst)/norm(oEst - pEst);
up = [0;0;1]; z0 = up - (up.'*x0)*x0; z0 = z0/norm(z0); y0 = cross(z0, x0);
R.T_E0_est = [x0 y0 z0 pEst; 0 0 0 1];
R.centre = oEst; R.normal = n; R.limbus_c = cEst; R.limbus_r = rEst; R.entry = pEst;
R.err_entry  = norm(pEst - P.eye.pE);
R.err_centre = norm(oEst);
R.err_axis   = acos(min(1, n(3)));
R.err_frame  = acos(min(1, max(-1, (trace(P.T_E0(1:3,1:3).'*R.T_E0_est(1:3,1:3)) - 1)/2)));
R.err_radius = rEst - rL;
R.X = Xs; R.Xtrue = Xt; R.uv = uv; R.cams = cams;
% ---------------- drawing data for the left image ----------------
th = linspace(0, 2*pi, 200);
R.img.limbus   = proj(cams(1).Pm, [rL*cos(th); rL*sin(th); zL + 0*th]);
R.img.pupil    = proj(cams(1).Pm, [1.92e-3*cos(th); 1.92e-3*sin(th); 10.26e-3 + 0*th]);
R.img.fit      = proj(cams(1).Pm, cEst + rEst*(e1*cos(th) + e2*sin(th)));
R.img.centre   = proj(cams(1).Pm, cEst);
R.img.entry    = proj(cams(1).Pm, pEst);
R.img.entryTrue = proj(cams(1).Pm, P.eye.pE);
d = cams(1).Cc; dist = norm(d); wv = d/dist;                % sphere silhouette seen by the left camera
a1 = cross(wv, [0;0;1]); a1 = a1/norm(a1); a2 = cross(wv, a1);
rs = P.eye.R*sqrt(1 - (P.eye.R/dist)^2); cs = (P.eye.R^2/dist)*wv;
R.img.globe = proj(cams(1).Pm, cs + rs*(a1*cos(th) + a2*sin(th)));
R.img.res = C.res;
R.fitCircle3D = cEst + rEst*(e1*cos(th) + e2*sin(th));
end

function uv = proj(Pm, X)
x = Pm*[X; ones(1, size(X,2))];
uv = x(1:2,:)./x(3,:);
end

function Rm = axang(w, a)
Kx = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
Rm = eye(3) + sin(a)*Kx + (1 - cos(a))*Kx*Kx;
end
