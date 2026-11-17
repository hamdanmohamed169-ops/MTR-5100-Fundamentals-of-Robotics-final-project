function H = fp_scene(ax, P, opts)
%FP_SCENE  Build the 3-D prototype scene: operating table, patient, eye and
%   the 3-DOF RCM injection robot on its 4-DOF arm, in axes ax. Units: mm, eye frame {E}.
%
%   Robot design (kinematically identical to the DH model used everywhere):
%     - stereo camera on a boom (locates the eye, FP_EYE_REGISTER)
%     - active 4-DOF positioning arm (closed-form IK, FP_ARM_IK) that carries
%       the whole RCM unit from its parking pose to the eye and back
%     - J1 yaw unit on z0 (axis passes through the RCM)
%     - yaw arc link centred on the RCM -> J2 pitch unit on z1
%     - pitch arc link centred on the RCM -> tool unit on the needle axis
%     - J3 insertion actuator + glass syringe + 30G needle (3-DOF R-R-P)
%   Each moving body is attached to an hgtransform that follows its DH frame
%   (body i -> T_E0*T0i, needle tip = origin of frame 3); see FP_SCENE_UPDATE.
%
%   opts.head  (true)  patient, table, drape, speculum
%   opts.setup (true)  4-DOF positioning arm (+ ghost arm for IK previews) + stereo camera
%   opts.xray  (false) transparent sclera (close-up view)
if nargin < 3, opts = struct(); end
if ~isfield(opts,'head'),  opts.head  = true;  end
if ~isfield(opts,'setup'), opts.setup = true;  end
if ~isfield(opts,'xray'),  opts.xray  = false; end
hold(ax,'on');
C.white  = [0.96 0.96 0.97]; C.silver = [0.74 0.76 0.79]; C.dark = [0.13 0.14 0.16];
C.rubber = [0.22 0.23 0.25]; C.steel  = [0.86 0.87 0.89]; C.led  = [0.20 0.75 1.00];
C.glass  = [0.85 0.93 1.00]; C.drape  = [0.25 0.46 0.58];

%% =========================== Eye ========================================
H.eye = hgtransform('Parent', ax);
R = P.eye.R*1e3; n = 96;
[sx,sy,sz] = sphere(n);
if opts.xray
    surf(R*sx, R*sy, R*sz, 'Parent', H.eye, 'FaceColor', [1 0.98 0.95], 'EdgeColor','none', ...
        'FaceAlpha', 0.22, 'FaceLighting','gouraud', 'AmbientStrength', 0.6);
else
    Xs = R*sx; Xs(R*sz > sqrt(R^2 - (P.eye.limbus_r*1e3)^2)) = NaN;     % opening for the cornea (iris visible)
    surf(Xs, R*sy, R*sz, 'Parent', H.eye, 'CData', scleraTexture(600, R), 'FaceColor','texturemap', ...
        'EdgeColor','none', 'FaceLighting','gouraud', 'AmbientStrength', 0.42, 'DiffuseStrength', 0.72, ...
        'SpecularStrength', 0.65, 'SpecularExponent', 55);
end
Rr = P.eye.Rret*1e3; Xr = Rr*sx; Xr(sz > 0.35) = NaN;
surf(Xr, Rr*sy, Rr*sz, 'Parent', H.eye, 'FaceColor', [0.85 0.30 0.26], 'EdgeColor','none', 'FaceAlpha', 0.30);
[dx,dy,dz] = sphere(20);
surf(0.9*dx, 0.9*dy, 0.15*dz - Rr+0.2, 'Parent', H.eye, 'FaceColor', [0.5 0.12 0.1], 'EdgeColor','none');
surf(0.9*dx + 3.5, 0.9*dy, 0.15*dz - Rr+0.6, 'Parent', H.eye, 'FaceColor', [1 0.9 0.6], 'EdgeColor','none');
L = P.eye.lens_ax*1e3; Lc = P.eye.lens_c*1e3;
surf(L(1)*sx+Lc(1), L(2)*sy+Lc(2), L(3)*sz+Lc(3), 'Parent', H.eye, 'FaceColor', [1 0.86 0.5], ...
    'EdgeColor','none', 'FaceAlpha', 0.45, 'SpecularStrength', 0.8);
% iris (procedural texture in polar coordinates) + pupil
[rr, tt] = ndgrid(linspace(0.33, 1, 72), linspace(0, 2*pi, 241));
ri = 5.75*rr; zi = 10.25 - 0.5*rr.^2 + 0.12*sin(pi*(rr - 0.33)/0.67);   % slightly domed iris
surf(ri.*cos(tt), ri.*sin(tt), zi, 'Parent', H.eye, 'CData', irisTexture(240, 720), ...
    'FaceColor', 'texturemap', 'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'AmbientStrength', 0.6, ...
    'SpecularStrength', 0.2);
th = linspace(0,2*pi,90);
[r1, t1] = ndgrid([0 1.95], th);                                         % pupil (black disc)
surf(r1.*cos(t1), r1.*sin(t1), 10.26 + 0*r1, 'Parent', H.eye, 'FaceColor', [0.01 0.01 0.01], ...
    'EdgeColor', 'none', 'SpecularStrength', 0, 'DiffuseStrength', 0, 'AmbientStrength', 0);   % stays black under any lighting
[r2, t2] = ndgrid([1.92 2.12], th);                                       % pupillary ruff
surf(r2.*cos(t2), r2.*sin(t2), 10.265 + 0*r2, 'Parent', H.eye, 'FaceColor', [0.18 0.10 0.05], 'EdgeColor', 'none');
if ~opts.xray                                                             % grey limbal ring
    [r3, t3] = ndgrid(linspace(5.55, 6.1, 4), th);
    surf(r3.*cos(t3), r3.*sin(t3), sqrt(R^2 - r3.^2) + 0.02, 'Parent', H.eye, 'FaceColor', [0.55 0.53 0.52], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.55, 'FaceLighting', 'gouraud');
end
% cornea (glossy, transparent)
Rc = 7.8; zc = 5.6; [cx,cy,cz] = sphere(70);
Zc = Rc*cz + zc; Xc = Rc*cx; Xc(Zc < 10.55) = NaN;
surf(Xc, Rc*cy, Zc, 'Parent', H.eye, 'FaceColor', [0.9 0.95 1], 'EdgeColor','none', 'FaceAlpha', 0.16, ...
    'SpecularStrength', 1, 'SpecularExponent', 120, 'DiffuseStrength', 0.08, 'AmbientStrength', 0.05);
if ~opts.xray
    catchlight(H.eye, [0; 0; zc], Rc, [-0.22; 0.30; 1], 0.95, 0.62, 0.92);   % operating-lamp reflection
    catchlight(H.eye, [0; 0; zc], Rc, [0.30; -0.22; 1], 0.30, 0.22, 0.70);   % secondary reflection
end
% keep-out zone + pars plana entry
K = P.eye.keepout_c*1e3; rk = P.eye.keepout_r*1e3;
surf(rk*dx+K(1), rk*dy+K(2), rk*dz+K(3), 'Parent', H.eye, 'FaceColor', [0.95 0.15 0.15], ...
    'EdgeColor','none', 'FaceAlpha', 0.85);
pE = P.eye.pE*1e3;
H.pEmark = plot3(pE(1), pE(2), pE(3), 'o', 'Color', [0.2 1 0.4], 'MarkerSize', 5, 'LineWidth', 1.5, 'Parent', H.eye);

%% =================== Patient, table, drape, speculum ===================
if opts.head
    [Fh, Vh, Nh, Ch, Bm] = fp_head_mesh();
    H.head = patch(ax, 'Faces', Fh, 'Vertices', Vh, 'VertexNormals', Nh, 'FaceVertexCData', Ch, ...
        'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.40, 'DiffuseStrength', 0.78, 'SpecularStrength', 0.10, 'SpecularExponent', 6, ...
        'BackFaceLighting', 'lit');
    % body under the surgical blanket (cloth folds from the same SDF method)
    H.body = patch(ax, 'Faces', Bm.F, 'Vertices', Bm.V, 'VertexNormals', Bm.N, 'FaceVertexCData', Bm.C, ...
        'FaceColor', 'interp', 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
        'AmbientStrength', 0.40, 'DiffuseStrength', 0.75, 'SpecularStrength', 0.04, 'BackFaceLighting', 'lit');
    % head rest (gel ring pillow)
    lathe(ax, [-32 70 -184], [-32 -70 -184], [0 0; 0 30; 4 34; 136 34; 140 30; 140 0], [0.30 0.55 0.62], 0.9);
    % operating table: padded mattress on a frame, chrome side rails with rail clamps
    rbox(ax, [-345 -440 -198], [175 210 -170], 10, [0.17 0.18 0.20]);
    rbox(ax, [-335 -432 -220], [165 202 -196], 4, [0.55 0.58 0.62]);
    for xr = [178 -348]
        tube(ax, [xr -430 -194], [xr 200 -194], 4, C.steel);
        for yr = [-400 -150 120], tube(ax, [xr yr -194], [xr - sign(xr)*12 yr -194], 3, C.steel); end
    end
    blob(ax, [-32 -40 -169.5], [80 150], 0.35);                         % contact shadow
    % table pedestal + base, and the operating-room floor
    zf = -560;
    rbox(ax, [-150 -260 -308], [-20 -120 -220], 10, [0.55 0.58 0.62]);
    lathe(ax, [-85 -190 -308], [-85 -190 zf+30], [0 0; 0 55; 222 50; 222 0], [0.86 0.87 0.89]);
    rbox(ax, [-300 -420 zf], [130 40 zf+30], 12, [0.55 0.58 0.62]);
    [fx, fy] = meshgrid(linspace(-1300, 1100, 2), linspace(-1300, 1100, 2));
    [gx, gy] = meshgrid(linspace(-1, 1, 64));
    fcd = 0.16*exp(-(gx.^2 + gy.^2)/0.35) + 0.07;
    surf(fx, fy, zf + 0*fx, 'Parent', ax, 'CData', repmat(fcd, 1, 1, 3).*reshape([0.80 0.90 1.0], 1, 1, 3), ...
        'FaceColor', 'texturemap', 'EdgeColor', 'none', 'FaceLighting', 'none');
    blob(ax, [-85 -190 zf+0.5], [380 380], 0.5);
    % Barraquer wire lid speculum: blades along both lid margins + spring loop temporally
    a = 13.0; b = 9.5; rv = (a^2/b + b)/2; cv = rv - b; Rl = 15.0;
    xb = linspace(-9.5, 9.5, 40);
    for s = [1 -1]
        yb = s*(-cv + sqrt(rv^2 - xb.^2) + 0.9);
        zb = sqrt(max(Rl^2 - xb.^2 - yb.^2, 0)) + 0.6;
        tubePath(ax, [xb; yb; zb], 0.45, C.steel);
        tubePath(ax, bezier([xb(end) yb(end) zb(end)], [17 s*6 zb(end)+1], [24 s*4.5 4], [30 s*3 3], 20), 0.45, C.steel);
    end
    tt2 = linspace(-pi/2, pi/2, 30);
    tubePath(ax, [30 + 4*cos(tt2); 3*sin(tt2); 3 + 0*tt2], 0.5, C.steel);
    % eyelashes: two staggered rows of tapered, curved lashes with random length, spacing and
    % angle (operated eye: pushed outwards by the speculum blades; fellow eye: closed, pointing down)
    rsL = RandStream('mt19937ar', 'Seed', 21);
    zl = @(xx, yy) sqrt(max(Rl^2 - xx.^2 - yy.^2, 0)) - 0.3;
    vy = @(xx) -cv + sqrt(rv^2 - xx.^2);
    xu = sort(-10.8 + 21.6*rand(rsL, 1, 66)); xd = sort(-10.2 + 20.4*rand(rsL, 1, 40));
    rowU = 0.25*(mod(1:numel(xu), 2) - 0.5); rowD = 0.25*(mod(1:numel(xd), 2) - 0.5);
    up = [xu; vy(xu) + 0.45 + rowU; zl(xu, vy(xu) + 0.45)];
    lo = [xd; -vy(xd) - 0.45 - rowD; zl(xd, vy(xd) + 0.45)];
    jit = @(n, s) s*randn(rsL, 1, n);
    [L1, A1] = lashSet(up, [0.13*xu/11 + jit(numel(xu), 0.08); 0.72 + jit(numel(xu), 0.06); 0.68 + 0*xu], ...
        (4.6 + 2.8*(1 - (xu/11).^2)).*(0.80 + 0.35*rand(rsL, 1, numel(xu))), [0; 0.30; 0.75], 0.075, rsL);
    [L2, A2] = lashSet(lo, [0.10*xd/11 + jit(numel(xd), 0.08); -0.80 + jit(numel(xd), 0.06); 0.48 + 0*xd], ...
        (2.8 + 1.8*(1 - (xd/10.5).^2)).*(0.80 + 0.35*rand(rsL, 1, numel(xd))), [0; -0.25; 0.30], 0.055, rsL);
    xf = sort(-74.5 + 21*rand(rsL, 1, 70)); yf = -2 + 0.012*(xf + 64).^2 + 0.2*(mod(1:numel(xf), 2) - 0.5);
    zf2 = sqrt(max(14.6^2 - (xf + 64).^2 - yf.^2, 0)) - 0.6;
    [L3, A3] = lashSet([xf; yf; zf2], [-0.13*(xf + 64)/11 + jit(numel(xf), 0.07); -0.86 + 0*xf; 0.44 + jit(numel(xf), 0.06)], ...
        (5.0 + 2.6*(1 - ((xf + 64)/11).^2)).*(0.80 + 0.35*rand(rsL, 1, numel(xf))), [0; -0.12; 0.42], 0.07, rsL);
    L = [L1; L2; L3]; A = [A1; A2; A3];
    H.lashes = patch(ax, 'Faces', reshape(1:size(L,1), 4, []).', 'Vertices', L, 'FaceColor', [0.06 0.045 0.035], ...
        'EdgeColor', 'none', 'FaceLighting', 'none', 'FaceVertexAlphaData', A, 'FaceAlpha', 'interp', ...
        'AlphaDataMapping', 'none');
    % tear meniscus along both lid margins (thin wet line where lid meets the globe)
    xm = linspace(-11.6, 11.6, 60);
    for s = [1 -1]
        ym = s*(vy(xm) - 0.15); zm = sqrt(max(12.15^2 - xm.^2 - ym.^2, 0));
        tubePath(ax, [xm; ym; zm], 0.16, [0.96 0.90 0.90], 0.55);
    end
    % caruncle + plica semilunaris at the medial canthus (towards the nose, -x)
    [ex, ey, ez] = sphere(18);
    surf(1.35*ex - 11.4, 1.05*ey, 0.75*ez + 5.2, 'Parent', ax, 'FaceColor', [0.86 0.52 0.52], 'EdgeColor', 'none', ...
        'FaceLighting', 'gouraud', 'SpecularStrength', 0.6, 'SpecularExponent', 30);
    tp = linspace(-1, 1, 20);
    tubePath(ax, [-10.4 + 0.3*tp.^2; 1.8*tp; sqrt(max(12.1^2 - (10.4 - 0.3*tp.^2).^2 - (1.8*tp).^2, 0)) + 0.05], 0.22, [0.90 0.70 0.68], 0.30);
end

%% =========================== Robot ======================================
TE0 = P.T_E0; TE0(1:3,4) = TE0(1:3,4)*1e3;
H.base = hgtransform('Parent', ax, 'Matrix', TE0);
H.L1 = hgtransform('Parent', ax); H.L2 = hgtransform('Parent', ax); H.L3 = hgtransform('Parent', ax);
H.arm = []; H.ghost = []; H.camHead = [];
H.plunger = hgtransform('Parent', H.L3);

% ---------- 4-DOF positioning arm (moves the whole RCM unit) + stereo camera ----------
if opts.setup
    H.arm   = buildArm(ax, P, 1, []);
    H.ghost = buildArm(ax, P, 0.18, [0.35 0.80 1.0]);        % IK-solution preview
    set(H.ghost.T, 'Visible', 'off');
    lathe(H.base, [0 0 172], [0 0 182], [0 0; 0 23; 2 24; 8 24; 10 23; 10 0], C.steel);   % quick-connect plate
    H.camHead = buildCamera(ax, P, C);
end

% ---------- J1 yaw unit (housing fixed in {0}) ----------
lathe(H.base, [0 0 104], [0 0 172], [0 0; 0 18; 2 20; 46 20; 54 18.5; 62 14; 66 8; 68 0], C.white);
lathe(H.base, [0 0 98],  [0 0 104], [0 0; 0 19.5; 6 19.5; 6 0], C.rubber);
lathe(H.base, [0 0 112], [0 0 115], [0 0; 0 20.4; 3 20.4; 3 0], C.led);
% ---------- body 1: yaw output flange + bent link + curved guide rail for J2 ----------
% J2 (pitch about z1) is realised as an arc-guide: a circular rail of radius Rg
% centred on the RCM, offset cg along z1 from the tool plane. The J2 carriage
% (body 2) rides on it, so no bearing sits on the pitch axis next to the face.
% Rail points in {1}: Rg*(-sin q2, cos q2, 0) + [0 0 cg],  q2 in [45, 135] deg.
Rg = 80; cg = 26; qg = linspace(deg2rad(42), deg2rad(138), 40);
rail = [-Rg*sin(qg); Rg*cos(qg); cg + 0*qg];
lathe(H.L1, [0 -98 0], [0 -92 0], [0 0; 0 18.5; 6 18.5; 6 0], C.dark);                     % J1 output flange
tubePath(H.L1, bezier([0 -92 0], [0 -78 4], [-46 -78 cg], rail(:,end).', 24), 9, C.white);   % bent link
sph(H.L1, rail(:,end), 9, C.white);
tubePath(H.L1, rail, 5, C.white);                                                          % guide rail
tubePath(H.L1, rail - [0; 0; 5.5], 1.6, C.steel);                                          % rack / encoder track
for e = [1 size(rail,2)], sph(H.L1, rail(:,e), 6.5, C.dark); end                            % end stops
lathe(H.L1, [0 -92 0], [0 -86 0], [0 0; 0 14; 6 12; 6 0], C.led);
% ---------- body 2: J2 carriage + tool unit + insertion actuator (pitch axis = +y2, needle = z2) ----------
rbox(H.L2, [-15 cg-9 -Rg-10], [15 cg+9 -Rg+10], 4, C.white);                                % carriage on the rail
lathe(H.L2, [0 cg+9 -Rg], [0 cg+30 -Rg], [0 0; 0 9; 18 9; 21 6; 21 0], C.white);              % J2 motor
lathe(H.L2, [0 cg+21 -Rg], [0 cg+23 -Rg], [0 0; 0 9.4; 2 9.4; 2 0], C.led);
rbox(H.L2, [-9 10 -Rg-8], [9 cg-8 -Rg+8], 3, C.silver);                                      % bracket to the tool unit
lathe(H.L2, [0 0 -58], [0 0 -104], [0 0; 0 13; 2 15.5; 40 15.5; 44 13; 46 0], C.white);         % tool unit
lathe(H.L2, [0 0 -54], [0 0 -58], [0 0; 0 15.5; 4 15.5; 4 0], C.rubber);
kn = [0 0; 0 7];
for k = 0:13, kn = [kn; 1+k 7+0.6*mod(k,2)]; end %#ok<AGROW>
kn = [kn; 15 7; 15 4; 18 3; 18 0];
lathe(H.L2, [0 0 -104], [0 0 -122], kn, C.steel);                                              % knurled cap
lathe(H.L2, [0 0 -122], [0 0 -128], [0 0; 0 3.6; 6 3.0; 6 0], C.dark);                          % cable strain relief
lathe(H.L2, [0 0 -54], [0 0 -46], [0 0; 0 13; 8 13; 8 0], C.steel);                            % nose collar
rbox(H.L2, [-6 -9 -54], [26 9 -47], 2, C.silver);                                              % clamp bracket
lathe(H.L2, [19 0 -56], [19 0 -18], [0 0; 0 5.5; 2 6.2; 36 6.2; 38 5.5; 38 0], C.steel);        % J3 insertion actuator
lathe(H.L2, [19 0 -18], [19 0 -14], [0 0; 0 3.5; 4 3.5; 4 0], C.dark);
lathe(H.L2, [19 0 -58], [19 0 -56], [0 0; 0 5; 2 5; 2 0], C.led);

% ---------- body 3: carriage + glass syringe + 30G needle (tip at origin) ----------
zh = -P.needle.L*1e3;
rbox(H.L3, [3 -4.5 zh-27], [17 4.5 zh-15], 1.5, C.silver);                  % carriage
lathe(H.L3, [0 0 zh-46], [0 0 zh-4], [0 0; 0 3.3; 42 3.3; 42 0], C.glass, 0.35);
for g = 0:6                                                                   % graduation marks
    lathe(H.L3, [0 0 zh-40+5*g], [0 0 zh-39.6+5*g], [0 0; 0 3.36; 0.4 3.36; 0.4 0], [0.25 0.3 0.35]);
end
lathe(H.L3, [0 0 zh-48], [0 0 zh-46], [0 0; 0 6; 2 6; 2 0], C.steel);       % finger flange
lathe(H.L3, [0 0 zh-4], [0 0 zh], [0 0; 0 1.6; 2 1.6; 4 0.6; 4 0], [0.15 0.45 0.95]);   % hub
lathe(H.L3, [0 0 zh], [0 0 0], [0 0; 0 0.16; abs(zh)-0.6 0.16; abs(zh) 0.02; abs(zh) 0], [0.92 0.92 0.95]);
lathe(H.plunger, [0 0 zh-40], [0 0 zh-38.5], [0 0; 0 3.1; 1.5 3.1; 1.5 0], [0.1 0.1 0.1]);
lathe(H.plunger, [0 0 zh-62], [0 0 zh-40], [0 0; 0 3.6; 2 3.6; 2 1.0; 22 1.0; 22 0], [0.96 0.96 0.96]);
patch('XData',[0 0.5 0], 'YData',[0 0 0], 'ZData',[0 -0.9 -0.9], 'FaceColor',[0.95 0.2 0.2], ...
    'EdgeColor','none','Parent',H.L3);                                       % bevel marker (set at mounting)
% ---------- power/signal cable: tool unit -> port on the J1 housing -> arm flange ----------
H.cablePort = [-21; 0; 140]; H.cableDir = [-1; 0; 0];                       % in {0}
lathe(H.base, [-20 0 140], [-26 0 140], [0 0; 0 4.2; 6 3.6; 6 0], C.steel);          % cable port
tubePath(H.base, bezier([-20 0 140], [-14 0 140], [-12 0 160], [-16 0 178], 14), 2.3, C.dark);   % into the quick-connect
H.tubeMesh = @tubeMesh; H.bez = @bezier;
H.cable = surf(ax, nan(2), nan(2), nan(2), 'FaceColor', C.dark, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
    'AmbientStrength', 0.38, 'DiffuseStrength', 0.72, 'SpecularStrength', 0.5, 'SpecularExponent', 15);
H.opts = opts;
fp_scene_update(H, P.q_home, P);
end

function [L, A] = lashSet(B, D, len, curl, w0, rs)
% eyelash ribbons: base points B (3xn), directions D (3xn), lengths len, curl vector,
% base half-width w0 (mm). Each lash: 6 segments, curl grows with t^2 (+ random amount),
% width tapers to a fine tip and the opacity fades towards the tip (per-vertex alpha A).
nS = 6; n = size(B,2); L = zeros(n*4*nS, 3); A = zeros(n*4*nS, 1); m = 0; t = linspace(0, 1, nS+1);
for i = 1:n
    d = D(:,i)/norm(D(:,i));
    cu = curl*(0.7 + 0.6*rand(rs));
    P = B(:,i) + len(i)*(d*t + cu*t.^2);
    sv = cross(d, [0;0;1]); if norm(sv) < 0.1, sv = [1;0;0]; end
    sv = sv/norm(sv);
    w = w0*(1 - 0.9*t.^0.8); al = 0.95*(1 - 0.8*t.^1.5);
    for k = 1:nS
        q = [P(:,k) - w(k)*sv, P(:,k) + w(k)*sv, P(:,k+1) + w(k+1)*sv, P(:,k+1) - w(k+1)*sv].';
        L(m+(1:4),:) = q; A(m+(1:4)) = [al(k) al(k) al(k+1) al(k+1)]; m = m + 4;
    end
end
end

function catchlight(parent, c0, Rc, dir, a, b, amax)
% soft elliptical specular highlight lying on the cornea sphere (centre c0, radius Rc)
dir = dir/norm(dir); p0 = c0 + Rc*dir;
e1 = cross([0;1;0], dir); e1 = e1/norm(e1); e2 = cross(dir, e1);
[rr, tt] = ndgrid(linspace(0, 1, 10), linspace(0, 2*pi, 36));
X = zeros(size(rr)); Y = X; Z = X;
for k = 1:numel(rr)
    q = p0 + a*rr(k)*cos(tt(k))*e1 + b*rr(k)*sin(tt(k))*e2;
    q = c0 + (Rc + 0.03)*(q - c0)/norm(q - c0);
    X(k) = q(1); Y(k) = q(2); Z(k) = q(3);
end
surf(X, Y, Z, 'Parent', parent, 'FaceColor', [1 1 1], 'EdgeColor', 'none', 'FaceLighting', 'none', ...
    'AlphaData', amax*(1 - rr.^2).^1.5, 'FaceAlpha', 'interp', 'AlphaDataMapping', 'none');
end

%% ============================ positioning arm ===========================
function A = buildArm(ax, P, alpha, tint)
% Each body is drawn in the arm-base frame at q = 0 (column up, the three
% planar links stretched out horizontally) under inner = inv(T_i(0))*T_Eb,
% so setting the outer transform to T_i(q) moves it exactly with its DH frame.
mmT = @(T) [T(1:3,1:3) T(1:3,4)*1e3; 0 0 0 1];
[~, T0all] = fp_arm_fk(zeros(4,1), P);
Teb = mmT(P.arm.T_Eb);
W = [0.95 0.95 0.96]; S = [0.74 0.76 0.79]; D = [0.14 0.15 0.17]; L = [0.20 0.75 1.00]; St = [0.86 0.87 0.89];
if ~isempty(tint), W = tint; S = tint; D = tint; L = tint; St = tint; end
A.T = gobjects(1,5); in = gobjects(1,5);
for i = 0:4
    A.T(i+1) = hgtransform('Parent', ax);
    if i == 0, Ti0 = eye(4); else, Ti0 = mmT(T0all(:,:,i)); end
    in(i+1) = hgtransform('Parent', A.T(i+1), 'Matrix', Ti0\Teb);
end
d1 = P.arm.d1*1e3; a2 = P.arm.a2*1e3; a3 = P.arm.a3*1e3; a4 = P.arm.a4*1e3;
xE = a2; xW = a2 + a3; xF = a2 + a3 + a4;
% body 0: rail clamp + column (static)
g = in(1);
rbox(g, [-34 -36 -22], [80 36 26], 8, D, alpha);
rbox(g, [52 -40 -30], [84 40 34], 6, S, alpha);
lathe(g, [0 0 -90], [0 0 0], [0 0; 0 40; 4 44; 86 44; 90 40; 90 0], W, alpha);
lathe(g, [0 0 0], [0 0 d1-75], [0 0; 0 44; 6 42; d1-90 38; d1-80 37; d1-75 34; d1-75 0], W, alpha);
lathe(g, [0 0 d1-84], [0 0 d1-80], [0 0; 0 38.5; 4 38.5; 4 0], L, alpha);
% body 1: turret (J1 base yaw) + shoulder housing (J2)
g = in(2);
lathe(g, [0 0 d1-75], [0 0 d1-58], [0 0; 0 37; 17 37; 17 0], D, alpha);
lathe(g, [0 0 d1-58], [0 0 d1-20], [0 0; 0 40; 30 40; 36 36; 38 0], W, alpha);
lathe(g, [0 -56 d1], [0 56 d1], filletCyl(112, 42, 8), W, alpha);
lathe(g, [0 -59 d1], [0 -56 d1], [0 0; 0 30; 3 30; 3 0], D, alpha);
lathe(g, [0 -60 d1], [0 -59 d1], [0 0; 0 22; 1 22; 1 0], L, alpha);
lathe(g, [0 56 d1], [0 59 d1], [0 0; 0 30; 3 30; 3 0], D, alpha);
% body 2: upper arm (a2) + elbow housing (J3)
g = in(3);
lathe(g, [0 0 d1], [xE 0 d1], [0 0; 0 30; 40 26; a2-40 24; a2 27; a2 0], W, alpha);
lathe(g, [60 0 d1], [64 0 d1], [0 0; 0 27; 4 27; 4 0], S, alpha);
lathe(g, [xE-64 0 d1], [xE-60 0 d1], [0 0; 0 25; 4 25; 4 0], S, alpha);
lathe(g, [xE -44 d1], [xE 44 d1], filletCyl(88, 33, 7), W, alpha);
lathe(g, [xE -47 d1], [xE -44 d1], [0 0; 0 24; 3 24; 3 0], D, alpha);
lathe(g, [xE -48 d1], [xE -47 d1], [0 0; 0 17; 1 17; 1 0], L, alpha);
% body 3: forearm (a3) + wrist housing (J4)
g = in(4);
lathe(g, [xE 0 d1], [xW 0 d1], [0 0; 0 25; 30 22; a3-30 20; a3 22; a3 0], W, alpha);
lathe(g, [xE+90 0 d1], [xE+94 0 d1], [0 0; 0 22.6; 4 22.6; 4 0], L, alpha);
lathe(g, [xW -34 d1], [xW 34 d1], filletCyl(68, 26, 6), W, alpha);
lathe(g, [xW -37 d1], [xW -34 d1], [0 0; 0 19; 3 19; 3 0], D, alpha);
lathe(g, [xW -38 d1], [xW -37 d1], [0 0; 0 13; 1 13; 1 0], L, alpha);
% body 4: wrist link (a4) + flange (tool interface)
g = in(5);
lathe(g, [xW 0 d1], [xF-12 0 d1], [0 0; 0 19; 20 17; a4-12 16; a4-12 0], W, alpha);
lathe(g, [xW+40 0 d1], [xW+43 0 d1], [0 0; 0 17.4; 3 17.4; 3 0], L, alpha);
lathe(g, [xF-12 0 d1], [xF 0 d1], [0 0; 0 18; 4 18; 6 24; 12 24; 12 0], St, alpha);
% ---- details: joint end caps, cable hoses, parting lines, E-stop, clamp knob ----
Y = [0.98 0.80 0.10]; Rd = [0.85 0.10 0.08];
if ~isempty(tint), Y = tint; Rd = tint; end
g = in(1);
lathe(g, [0 -37 250], [0 -48 250], [0 0; 0 13; 11 13; 11 0], Y, alpha);                    % emergency stop
lathe(g, [0 -48 250], [0 -56 250], [0 0; 0 9; 3 9; 6 7; 8 4; 8 0], Rd, alpha);
tube(g, [84 0 0], [98 0 0], 4, St, alpha); sph(g, [104 0 0], 9, D, alpha);                  % rail-clamp knob
rbox(g, [54 -32 -84], [82 32 -24], 6, S, alpha);                                           % clamp jaw around the table rail
tube(g, [82 0 -60], [96 0 -60], 3.5, St, alpha); sph(g, [100 0 -60], 7, D, alpha);
g = in(2);
capRing(g, [0 -60 d1], [0 -1 0], 22, alpha, tint); capRing(g, [0 59 d1], [0 1 0], 30, alpha, tint);
g = in(3);
lathe(g, [xE 44 d1], [xE 47 d1], [0 0; 0 24; 3 24; 3 0], D, alpha);
capRing(g, [xE -48 d1], [0 -1 0], 17, alpha, tint); capRing(g, [xE 47 d1], [0 1 0], 24, alpha, tint);
lathe(g, [a2/2-1 0 d1], [a2/2+1 0 d1], [0 0; 0 25.4; 2 25.4; 2 0], D, alpha);               % parting line
tubePath(g, bezier([30 22 d1+20], [100 33 d1+27], [a2-100 33 d1+25], [a2-35 22 d1+18], 30), 4, D, alpha);
g = in(4);
lathe(g, [xW 34 d1], [xW 37 d1], [0 0; 0 19; 3 19; 3 0], D, alpha);
capRing(g, [xW -38 d1], [0 -1 0], 13, alpha, tint); capRing(g, [xW 37 d1], [0 1 0], 19, alpha, tint);
lathe(g, [xE+a3/2-1 0 d1], [xE+a3/2+1 0 d1], [0 0; 0 21.4; 2 21.4; 2 0], D, alpha);
tubePath(g, bezier([xE+30 18 d1+18], [xE+90 27 d1+22], [xW-90 26 d1+20], [xW-28 17 d1+15], 30), 3.4, D, alpha);
end

function Hc = buildCamera(ax, P, C)
% Stereo tracking camera (static) on a mobile floor stand behind the patient's head:
% weighted 5-leg base with casters, column, two-segment articulated boom with
% joint knuckles, gimbal yoke, and a camera head with two lenses (baseline
% P.cam.base, optical centres exactly at the model's camera centres) each
% surrounded by an LED ring light.
mm = 1e3; K = P.cam;
Cc = K.C*mm; tg = K.target*mm;
zc = (tg - Cc)/norm(tg - Cc);
bx = cross(zc, [0;0;1]); bx = bx/norm(bx); by = cross(zc, bx);
pole = [-300; 300]; zf = -560;                                     % floor level under the table
grey = [0.80 0.82 0.85];
% mobile base: hub, 5 legs, casters
lathe(ax, [pole; zf+40], [pole; zf+95], [0 0; 0 60; 10 58; 45 40; 55 30; 55 0], C.white);
for k = 0:4
    a = 2*pi*k/5 + 0.3; e = [cos(a); sin(a); 0];
    tip = [pole; zf+40] + 230*e;
    tubePath(ax, bezier([pole; zf+70], [pole; zf+70] + 80*e, tip + [0;0;20], tip, 12), 11, grey);
    lathe(ax, tip, tip - [0;0;14], [0 0; 0 12; 14 12; 14 0], C.dark);
    lathe(ax, tip - [0;0;26] - 10*[sin(a); -cos(a); 0], tip - [0;0;26] + 10*[sin(a); -cos(a); 0], ...
        [0 0; 0 13; 20 13; 20 0], C.rubber);                             % caster wheel
end
% column with a height-adjust collar and status ring
lathe(ax, [pole; zf+95], [pole; 640], [0 0; 0 30; 8 27; 640-zf-103 25; 640-zf-95 20; 640-zf-95 0], C.white);
lathe(ax, [pole; 120], [pole; 170], filletCyl(50, 33, 6), grey);
lathe(ax, [pole; 548], [pole; 556], [0 0; 0 27.6; 8 27.6; 8 0], C.led);
% articulated boom: shoulder knuckle -> segment 1 -> elbow knuckle -> segment 2 -> gimbal
J0 = [pole; 640]; M = Cc - 125*zc; J1 = [M(1:2) + 0.45*(pole - M(1:2)); 700];
knuckle = @(p) lathe(ax, p - [0;0;28], p + [0;0;28], filletCyl(56, 30, 8), C.white);
knuckle(J0); knuckle(J1);
capsuleLink(ax, J0 + [0;0;10], J1 + [0;0;10], 19, C.white);
capsuleLink(ax, J1 - [0;0;10], M + [0;0;30], 16, C.white);
lathe(ax, J0 + [0;0;28], J0 + [0;0;31], [0 0; 0 22; 3 22; 3 0], C.led);
sph(ax, M + [0;0;30], 20, grey);                                   % gimbal
% camera head (in its own frame: x = baseline, z = optical axis toward the eye)
Hc = hgtransform('Parent', ax, 'Matrix', [bx by zc Cc; 0 0 0 1]);
rbox(Hc, [-98 -34 -62], [98 34 -4], 12, C.white);                 % housing
rbox(Hc, [-96 -32 -8], [96 32 2], 6, C.dark);                     % front bezel
rbox(Hc, [-22 -20 -75], [22 20 -60], 5, grey);                    % mount block
for f = -3:3                                                        % cooling fins on the back
    rbox(Hc, [12*f-2 -26 -66], [12*f+2 26 -61], 1, grey);
end
tubePath(Hc, [[0;0;-74] (Hc.Matrix(1:3,1:3).'*(M + [0;0;10] - Cc))], 9, grey);   % yoke to the gimbal
b2 = K.base/2*mm;
for xL = [-b2 b2]
    c = [xL; 0; 0];
    lathe(Hc, c + [0;0;-2], c + [0;0;20], [0 0; 0 17; 16 17; 18 14; 22 14; 22 0], [0.08 0.08 0.09]);   % lens barrel
    lathe(Hc, c + [0;0;6], c + [0;0;8], [0 0; 0 17.5; 2 17.5; 2 0], [0.45 0.46 0.48]);                % focus ring
    lathe(Hc, c + [0;0;19.5], c + [0;0;20.5], [0 0; 0 10; 1 10; 1 0], [0.12 0.22 0.45]);              % front glass
    sph(Hc, c + [-3; -3; 20.6], 1.6, [0.9 0.95 1], 0.8);                                              % glass highlight
    th = linspace(0, 2*pi, 13); th(end) = [];
    for t = th, sph(Hc, c + [24*cos(t); 24*sin(t); 2.5], 2.6, [0.92 0.97 1]); end                     % LED ring light
end
sph(Hc, [0; 24; 2.5], 2.2, C.led);                                 % status LED
end

function capsuleLink(parent, a, b, r, col)
a = a(:); b = b(:); Lb = norm(b - a);
lathe(parent, a, b, [0 0; 0 r; Lb r; Lb 0], col); sph(parent, a, r, col); sph(parent, b, r, col);
end

%% ============================ textures ==================================
function T = irisTexture(nr, nt)
% procedural brown iris (rows = radius pupil->limbus, cols = angle):
% radial stroma fibres, golden pupillary zone, collarette, crypts,
% circular contraction furrows and a dark limbal ring
rng(11, 'twister');
r = linspace(0, 1, nr).'; t = linspace(0, 1, nt);
fib = zeros(nr, nt);
for k = 1:9
    fib = fib + randn*sin(2*pi*(12*k + randi(30))*t + 2*pi*rand + 4*r*rand) / sqrt(k);
end
fib = 0.5 + 0.15*fib;
base = [0.34 0.21 0.10]; gold = [0.68 0.50 0.24]; dark = [0.09 0.055 0.03];
pz = exp(-((r - 0.08)/0.16).^2);                  % pupillary zone (golden)
rcol = 0.30 + 0.025*sin(2*pi*14*t + 1.3) + 0.012*sin(2*pi*37*t);   % zig-zag collarette radius
col = exp(-((r - rcol)/0.045).^2);                % collarette (raised, lighter)
fl = zeros(nr, nt);                               % light flecks in the stroma
for k = 1:60
    rc = 0.35 + 0.5*rand; tc = rand; dt_ = mod(t - tc + 0.5, 1) - 0.5;
    fl = fl + exp(-((r - rc)/0.012).^2 - (dt_/0.004).^2);
end
fib = fib + 0.25*min(fl, 1);
fur = 0.5*exp(-((r - 0.62)/0.015).^2) + 0.4*exp(-((r - 0.78)/0.015).^2);   % contraction furrows
lim = 1./(1 + exp(-(r - 0.9)/0.025));             % limbal ring
cry = zeros(nr, nt);                              % crypts: small dark lacunae
for k = 1:45
    rc = 0.3 + 0.5*rand; tc = rand; s1 = 0.02 + 0.03*rand; s2 = 0.006 + 0.01*rand;
    dt_ = mod(t - tc + 0.5, 1) - 0.5;
    cry = cry + exp(-((r - rc)/s1).^2 - (dt_/s2).^2);
end
cry = min(cry, 1);
T = zeros(nr, nt, 3);
for c = 1:3
    v = base(c) + (gold(c) - base(c))*min(1, 0.55*pz + 0.6*col) + 0.35*(fib - 0.5);
    v = v.*(1 - 0.55*cry).*(1 - 0.35*fur);
    T(:,:,c) = v.*(1 - lim) + dark(c)*lim;
end
T = min(max(T, 0), 1);
end

function T = scleraTexture(n, Rg)
% sclera + bulbar conjunctiva (rows = polar angle -z..+z, cols = azimuth):
% warm white near the limbus, pinker toward the fornix, branching episcleral
% vessels running toward the limbus
rng(5, 'twister');
row = linspace(0, 1, n).';
pink = 1./(1 + exp((row - 0.74)/0.03));           % 0 at the limbus side, 1 toward the equator
base = [0.96 0.94 0.91]; fornix = [0.93 0.80 0.77];
T = zeros(n, n, 3);
for c = 1:3, T(:,:,c) = repmat(base(c) + (fornix(c) - base(c))*pink, 1, n); end
T = T .* (1 + 0.012*randn(n, n));
for v = 1:60
    walk(round(n*(0.64 + 0.08*rand)), randi(n), 0.16 + 0.22*rand, 2);
end
% soft shadow of the lid margins on the globe: texel (i,j) of sphere(n-1) lies at
% phi = -pi/2 + pi*(i-1)/(n-1), theta = -pi + 2*pi*(j-1)/(n-1)
if nargin > 1
    [ph, tq] = ndgrid(linspace(-pi/2, pi/2, n), linspace(-pi, pi, n));
    gx = Rg*cos(ph).*cos(tq); gy = Rg*cos(ph).*sin(tq);
    a = 13.0; b = 9.5; rv = (a^2/b + b)/2; cv = rv - b;
    ap = max(sqrt(gx.^2 + (gy + cv).^2) - rv, sqrt(gx.^2 + (gy - cv).^2) - rv);   % vesica SDF (fissure)
    sh = 1 - 0.30*exp(min(ap, 0)/1.4) .* (ph > 0);
    T = T .* sh;
    T(:,:,1) = T(:,:,1) + 0.03*exp(min(ap, 0)/2.5) .* (ph > 0);                   % slightly pinker near the lids
end
T = min(max(T, 0), 1);
    function walk(r, c, strength, depth)
        len = round(n*(0.06 + 0.10*rand)); w = strength;
        for s = 1:len
            if r < 1 || r > n || row(r) > 0.835, break; end
            cc = mod(round(c) - 1, n) + 1;
            a = w*(0.55 + 0.45*pink(r));
            T(r, cc, :) = (1 - a)*T(r, cc, :) + a*reshape([0.80 0.36 0.34], 1, 1, 3);
            if w > 0.7 && cc < n
                T(r, cc+1, :) = (1 - a/2)*T(r, cc+1, :) + a/2*reshape([0.82 0.40 0.38], 1, 1, 3);
            end
            r = r + 1; c = c + randn*0.9; w = w*0.985;
            if depth > 0 && rand < 0.03, walk(r, c, 0.7*w, depth - 1); end
        end
    end
end

%% ============================ geometry helpers ==========================
function lathe(parent, a, b, prof, col, alpha)
% surface of revolution about the segment a->b; prof = [s r] (mm from a)
if nargin < 6, alpha = 1; end
a = a(:); b = b(:); d = b - a;
Rm = frameFromZ(d/norm(d));
th = linspace(0, 2*pi, 56);
X = prof(:,2)*cos(th); Y = prof(:,2)*sin(th); Z = repmat(prof(:,1), 1, numel(th));
pts = Rm*[X(:) Y(:) Z(:)].' + a;
surf(reshape(pts(1,:), size(X)), reshape(pts(2,:), size(X)), reshape(pts(3,:), size(X)), ...
    'Parent', parent, 'FaceColor', col, 'EdgeColor', 'none', 'FaceAlpha', alpha, ...
    'FaceLighting', 'gouraud', 'AmbientStrength', 0.38, 'DiffuseStrength', 0.72, ...
    'SpecularStrength', 0.85, 'SpecularExponent', 22, 'BackFaceLighting', 'lit');
end

function p = filletCyl(len, r, f)
% profile of a cylinder of length len, radius r, fillet radius f at both ends
a = linspace(0, pi/2, 7).';
p1 = [f - f*cos(a), r - f + f*sin(a)];
p2 = [len - f + f*sin(a), r - f + f*cos(a)];
p = [0 0; p1; p2; len 0];
end

function tube(parent, a, b, r, col, alpha)
if nargin < 6, alpha = 1; end
a = a(:); b = b(:); Lb = norm(b - a);
lathe(parent, a, b, [0 0; 0 r; Lb r; Lb 0], col, alpha);
end

function sph(parent, c, r, col, alpha)
if nargin < 5, alpha = 1; end
[X,Y,Z] = sphere(36); c = c(:);
surf(r*X + c(1), r*Y + c(2), r*Z + c(3), 'Parent', parent, 'FaceColor', col, 'EdgeColor','none', 'FaceAlpha', alpha, ...
    'FaceLighting','gouraud', 'AmbientStrength', 0.38, 'DiffuseStrength', 0.72, 'SpecularStrength', 0.85, 'SpecularExponent', 22);
end

function capsule(parent, a, b, r, col)
a = a(:); b = b(:); Lb = norm(b - a);
lathe(parent, a, b, [0 0; 0 r; Lb r; Lb 0], col); sph(parent, a, r, col); sph(parent, b, r, col);
u = (b - a)/Lb;                                     % thin silver accent bands
lathe(parent, a + u*r*1.3, a + u*(r*1.3+2), [0 0; 0 r+0.4; 2 r+0.4; 2 0], [0.74 0.76 0.79]);
lathe(parent, b - u*(r*1.3+2), b - u*r*1.3, [0 0; 0 r+0.4; 2 r+0.4; 2 0], [0.74 0.76 0.79]);
end

function arcTube(parent, rad, e1, e2, a0, a1, r, col, n)
ang = linspace(a0, a1, n+1);
pts = rad*(e1*cos(ang) + e2*sin(ang));
tubePath(parent, pts, r, col);
sph(parent, pts(:,1), r, col); sph(parent, pts(:,end), r, col);
end

function tubePath(parent, pts, r, col, alpha)
% smooth tube (one surface) along a 3-D polyline pts (3 x n), parallel transport frame
if nargin < 5, alpha = 1; end
[X, Y, Z] = tubeMesh(pts, r);
surf(X, Y, Z, 'Parent', parent, 'FaceColor', col, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', ...
    'FaceAlpha', alpha, 'AmbientStrength', 0.38, 'DiffuseStrength', 0.72, 'SpecularStrength', 0.85, 'SpecularExponent', 22);
end

function [X, Y, Z] = tubeMesh(pts, r)
n = size(pts, 2); th = linspace(0, 2*pi, 28);
X = zeros(n, numel(th)); Y = X; Z = X; nPrev = [];
for k = 1:n
    if k == 1, t = pts(:,2) - pts(:,1); elseif k == n, t = pts(:,n) - pts(:,n-1);
    else, t = pts(:,k+1) - pts(:,k-1); end
    t = t/norm(t);
    if isempty(nPrev)
        nv = cross(t, [0;0;1]); if norm(nv) < 0.1, nv = cross(t, [0;1;0]); end
    else
        nv = nPrev - (nPrev.'*t)*t;
    end
    nv = nv/norm(nv); bv = cross(t, nv);
    ring = pts(:,k) + r*(nv*cos(th) + bv*sin(th));
    X(k,:) = ring(1,:); Y(k,:) = ring(2,:); Z(k,:) = ring(3,:);
    nPrev = nv;
end
end

function capRing(parent, c, ax, r, alpha, tint)
% joint end cap: brushed disc and 6 screw heads on a bolt circle
ax = ax(:)/norm(ax); c = c(:); Rm = frameFromZ(ax);
c1 = [0.80 0.82 0.85]; c2 = [0.35 0.36 0.38];
if ~isempty(tint), c1 = tint; c2 = tint; end
lathe(parent, c, c + 1.2*ax, [0 0; 0 0.62*r; 1.2 0.55*r; 1.2 0], c1, alpha);
for k = 0:5
    p = c + Rm*[0.78*r*cos(k*pi/3 + pi/6); 0.78*r*sin(k*pi/3 + pi/6); 0.3];
    sph(parent, p, 0.07*r, c2, alpha);
end
end

function pts = bezier(p0, p1, p2, p3, n)
t = linspace(0, 1, n);
pts = p0(:)*(1-t).^3 + 3*p1(:)*((1-t).^2.*t) + 3*p2(:)*((1-t).*t.^2) + p3(:)*t.^3;
end

function rbox(parent, lo, hi, rr, col, alpha)
if nargin < 6, alpha = 1; end
% rounded box (superellipsoid)
c = (lo(:) + hi(:))/2; s = (hi(:) - lo(:))/2;
e = max(0.08, min(0.45, rr/min(s)));
[u, v] = meshgrid(linspace(-pi/2, pi/2, 24), linspace(-pi, pi, 48));
f = @(w, ee) sign(w).*abs(w).^ee;
X = c(1) + s(1)*f(cos(u), e).*f(cos(v), e);
Y = c(2) + s(2)*f(cos(u), e).*f(sin(v), e);
Z = c(3) + s(3)*f(sin(u), e);
surf(X, Y, Z, 'Parent', parent, 'FaceColor', col, 'EdgeColor', 'none', 'FaceLighting', 'gouraud', 'FaceAlpha', alpha, ...
    'AmbientStrength', 0.38, 'DiffuseStrength', 0.7, 'SpecularStrength', 0.6, 'SpecularExponent', 15);
end

function blob(parent, c, ab, alpha)
% soft contact shadow (concentric translucent ellipses)
th = linspace(0, 2*pi, 60);
for k = 1:4
    s = 1 - 0.18*(k-1);
    patch('XData', c(1) + ab(1)*s*cos(th), 'YData', c(2) + ab(2)*s*sin(th), ...
        'ZData', c(3)*ones(size(th)) + 0.05*k, 'FaceColor', [0 0 0], 'FaceAlpha', alpha/4, ...
        'EdgeColor', 'none', 'Parent', parent, 'FaceLighting', 'none');
end
end

function Rm = frameFromZ(z)
tmp = [1;0;0]; if abs(z.'*tmp) > 0.9, tmp = [0;1;0]; end
x = cross(tmp, z); x = x/norm(x); y = cross(z, x); Rm = [x y z];
end
