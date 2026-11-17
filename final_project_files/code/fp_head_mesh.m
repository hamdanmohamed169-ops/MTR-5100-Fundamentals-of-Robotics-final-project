function [F, V, N, C, B, S] = fp_head_mesh()
%FP_HEAD_MESH  Patient (supine) for the 3-D prototype scene.
%   Head = signed-distance model triangulated with isosurface:
%     - the FRONT of the face is a height surface z = h(u,v) built from an
%       anthropometric mid-line profile (glabella, nasion, nose, lips, chin)
%       plus a lateral fall-off and smooth features (nose dorsum, tip and
%       alae, orbits, brow ridge, cheekbones, philtrum, vermilion with
%       cupid's bow, commissures, nasolabial folds, labiomental fold, chin)
%     - the face outline (front view) and the jaw underside close the face
%       volume, which is blended with the cranium, ears and neck
%     - eyelids: operated eye held open by the speculum, fellow eye closed
%     - surgical cap; the body under the blanket is a second, coarser SDF mesh.
%   Units mm, eye frame {E}: operated right eye centred at the origin,
%   zE up (out of the face), xE temporal, yE toward the forehead.
%   Proportions: adult anthropometry (face height nasion-menton ~125 mm,
%   bizygomatic width ~140 mm, IPD 64 mm, nasal tip ~28 mm in front of the
%   cornea, lips 5-7 mm behind the nose-chin line).
%   Outputs: faces F, vertices V, vertex normals N, vertex colours C (head),
%            B.F, B.V, B.N, B.C (body + blanket),
%            S.x, S.y, S.z (grid vectors), S.d  head SDF (meshgrid layout, used
%            for collision checks; normalised by its gradient so it never
%            over-estimates the distance), S.body = same for body + blanket.
%   Cached per MATLAB session.
persistent cache
if ~isempty(cache)
    F = cache.F; V = cache.V; N = cache.N; C = cache.C; B = cache.B; S = cache.S; return
end
%% ------------------------------- head ----------------------------------
h = 1.0;
[x, y, z] = meshgrid(single(-122:h:62), single(-240:h:136), single(-172:h:52));
mx = -32;                                                         % facial midline
d = headSdf(x, y, z, h, mx);
S = struct('x', squeeze(x(1,:,1)), 'y', squeeze(y(:,1,1)).', 'z', squeeze(z(1,1,:)).', 'd', d);   % grid vectors + SDF (meshgrid layout)
[F, V] = isosurface(x, y, z, d, 0);
V = double(V);
[F, V] = largestPart(F, V);                                       % drop tiny stray fragments
N = double(isonormals(x, y, z, d, V));
N = -N ./ max(vecnorm(N, 2, 2), eps);                             % outward normals
% Taubin smoothing (no shrinkage): light, a little stronger at the thin lid margins
w = 0.30 + 0.60*exp(-(sqrt(sum(V.^2, 2))/20).^4);
V = taubin(F, V, w, 8);
N2 = vertexNormal(triangulation(F, V));
if mean(sum(N2.*N, 2)) < 0, N2 = -N2; end
N = N2;
ao = ambientOcc(x, y, z, d, V, N, [1.5 4 9]);
% ---- vertex colours ----
uu = V(:,1) - mx; vv = V(:,2); front = 1./(1 + exp(-(V(:,3) + 5)/3));   % 1 on the face, 0 on the back of the head
skin = [0.86 0.67 0.57];
C = repmat(skin, size(V,1), 1);
gb = @(c, s) exp(-sum(((V - c)./s).^2, 2));                       % gaussian blob helper
hsh = @(a, b) sin(V(:,1)*a(1) + V(:,2)*a(2) + V(:,3)*a(3) + b);
mot = hsh([1.31 0.73 0.2], 0) .* hsh([0.4 1.17 0.91], 1) .* hsh([0.37 0.2 1.23], 2);
C = C .* (1 + 0.020*mot + 0.018*sin(V(:,1)/9 + 1.3*sin(V(:,2)/13)));     % mottling + broad tone
C = C + (gb([mx+38 -28 8], [15 14 14]) + gb([mx-38 -28 8], [15 14 14])) .* [0.05 -0.02 -0.025];   % cheeks
C = C + gb([mx -38 40], [8 7 9]) .* [0.05 -0.015 -0.025];                                        % nose tip
C = C + (gb([mx+13 -42 24], [6 6 6]) + gb([mx-13 -42 24], [6 6 6])) .* [0.05 -0.02 -0.03];       % alae
for s = [-1 1]
    C = C + gb([mx+s*76 -12 -52], [10 32 22]) .* [0.05 -0.02 -0.03];                            % ears
    C = C + gb([mx+s*30 -11 9], [12 4 6]) .* [-0.04 -0.05 -0.02];                               % under-eye shadow
    C = C + gb([mx+s*32 7 10], [12 4 5]) .* [0.01 -0.03 -0.03];                                 % upper-lid tone
end
C = C + gb([mx 45 12], [45 20 15]) .* [0.02 0.02 -0.01];                                         % forehead
% shaved beard shadow with stubble speckle
[lipU, lipL] = lipMasks(uu, vv);
lipM = max(lipU, lipL) .* front;
beard = (1./(1 + exp((vv + 50)/3)) .* exp(-(uu/50).^6) + 1./(1 + exp((vv + 92)/5)) .* exp(-(uu/62).^4)) ...
        .* (1 - exp(-((vv + 48)/4).^2) .* exp(-(uu/12).^2)) .* front;
beard = min(beard, 1) .* (1 - lipM);
stub = 0.75 + 0.25*hsh([7.1 5.3 3.7], 0.3).*hsh([4.3 6.7 5.9], 1.1);
C = C + (beard .* stub) .* [-0.045 -0.035 -0.01];
lipC = [0.72 0.43 0.42] + lipU.*[-0.03 -0.02 -0.01];
C = C.*(1 - lipM) + lipC.*lipM;                                                                  % lips
C = C + front .* exp(-((vv + 67)/0.9).^2) .* max(1 - (uu/25).^2, 0) .* [-0.28 -0.20 -0.17];      % mouth line
% eyebrows: arched band, dense at the head, thin at the tail
for s = [-1 1]
    t = (s*uu - 10)/44;
    vc = 15.5 + 7*sin(pi*min(max(t, 0), 1)*0.85) - 1.5*t;
    th = 4.4 - 2.4*min(max(t, 0), 1);
    br = 1./(1 + exp(-t/0.06)) .* 1./(1 + exp((t - 1)/0.08)) .* exp(-((vv - vc)./th).^4) .* front;
    br = br .* (0.90 + 0.10*hsh([0.9 0.3 0.5], s));
    C = C .* (1 - 0.68*br) + [0.24 0.17 0.12] .* (0.68*br);
end
% closed fellow eye: slightly pinker lid, dark lash line
xl = (V(:,1) + 64)/13;
lidL = exp(-abs(xl).^6) .* (V(:,3) > 0);
C = C + lidL .* exp(-((V(:,2) - 3)/5).^2) .* [0.02 -0.02 -0.01];
ll = 0.85*lidL.*exp(-((V(:,2) + 2 - 0.012*(V(:,1)+64).^2)/0.9).^2);
C = C .* (1 - ll) + ll .* [0.10 0.07 0.06];
% operated eye: pink conjunctiva / lid margins inside the fissure, slight redness around
rv = sqrt(sum(V.^2, 2));
conj = 1./(1 + exp((rv - 13.3)/0.3));
C = C .* (1 - conj) + [0.86 0.52 0.50] .* conj;
C = C + exp(-((rv - 15.5)/3).^2) .* (V(:,3) > 0) .* [0.03 -0.03 -0.03];
C = C + exp(-((V(:,2) - 11.5)/1.0).^2) .* (abs(V(:,1)) < 14) .* (V(:,3) > 4) .* [-0.06 -0.06 -0.04];
% surgical cap
capC = [0.30 0.50 0.62] .* (1 + 0.05*sin(V(:,1)/3 + V(:,3)/4).*sin(V(:,2)/5));
capV = 1./(1 + exp(-(V(:,2) - 66)/1.2));
C = C.*(1 - capV) + capC.*capV;
nV = size(V,1);                                                   % light colour smoothing over the mesh
Ah = sparse(F(:,[1 2 3 2 3 1]), F(:,[2 3 1 1 2 3]), 1, nV, nV); Ah = double(Ah > 0) + 2*speye(nV);
for it = 1:2, C = (Ah*C)./full(sum(Ah, 2)); end
C = C .* (0.60 + 0.40*ao);
C = min(max(C, 0), 1);
clear x y z d

%% ------------------------- body + blanket ------------------------------
% defined in the original table coordinates zo = z + 48 (head rest and
% mattress lowered 48 mm so the head has its real depth)
hb = 3.5;
[x, y, z] = meshgrid(single(-250:hb:200), single(-446:hb:-180), single(-172:hb:-28));
zo = z + 48;
db = sdEll(x,y,zo, [mx -340 -82], [150 150 52]);                         % chest / abdomen
db = smin(db, sdCaps(x,y,zo, [mx-130 -250 -84], [mx+130 -250 -84], 46), 40);   % shoulders
for s = [-1 1]
    db = smin(db, sdCaps(x,y,zo, [mx+s*165 -262 -90], [mx+s*172 -450 -96], 33), 30);   % arms
end
db = smin(db, sdCaps(x,y,zo, [mx -200 -62], [mx -260 -78], 60), 25);     % neck root / collar
fold = 2.4*sin(x/21 + 0.6*sin(y/37)) .* sin(y/29 + x/90) + 1.2*sin(x/9 + y/13);
dk = db - 9 + fold;
dk = smin(dk, max(max(abs(x - mx) - 200, zo + 100), -(zo + 121)), 30); % skirt falling onto the table
dk = smax(dk, -(zo + 121), 2);
dk = smax(dk, y + 205 - 0.25*(zo + 40), 6);                             % hem below the neck (slanted)
dk = smax(dk, -(y + 438), 3);                                           % foot end of the table
shell = max(dk, -(dk + 6));                                             % 6 mm thick blanket
db = smax(db, -(zo + 121), 2); db = smax(db, -(y + 438), 3);
da = smin(db, shell, 3);                                                % gown-covered body + blanket
[B.F, B.V] = isosurface(x, y, z, da, 0);
B.V = double(B.V);
[B.F, B.V] = largestPart(B.F, B.V);
N0 = double(isonormals(x, y, z, da, B.V)); N0 = -N0 ./ max(vecnorm(N0, 2, 2), eps);
B.V = taubin(B.F, B.V, ones(size(B.V,1), 1), 6);
B.N = vertexNormal(triangulation(B.F, B.V));
if mean(sum(B.N.*N0, 2)) < 0, B.N = -B.N; end
aob = ambientOcc(x, y, z, da, B.V, B.N, [4 10 20]);
B.C = repmat([0.22 0.47 0.55], size(B.V,1), 1) .* (1 + 0.03*sin(B.V(:,1)/4 + B.V(:,2)/6));
B.C = B.C + exp(-((B.V(:,2) + 203)/5).^2) .* [0.06 0.07 0.07];          % folded hem
vq = num2cell(single(B.V), 1);
gown = abs(double(interp3(x, y, z, db, vq{:}))) < abs(double(interp3(x, y, z, shell, vq{:})));
B.C(gown,:) = repmat([0.66 0.74 0.82], nnz(gown), 1);                   % hospital gown
nB = size(B.V,1);                                                      % smooth the gown / blanket colour boundary
Ab = sparse(B.F(:,[1 2 3 2 3 1]), B.F(:,[2 3 1 1 2 3]), 1, nB, nB); Ab = double(Ab > 0) + speye(nB);
for it = 1:4, B.C = (Ab*B.C)./full(sum(Ab, 2)); end
B.C = min(max(B.C .* (0.5 + 0.5*aob), 0), 1);
[g1, g2, g3] = gradient(da, hb); da = da ./ max(sqrt(g1.^2 + g2.^2 + g3.^2), 1);
S.body = struct('x', squeeze(x(1,:,1)), 'y', squeeze(y(:,1,1)).', 'z', squeeze(z(1,1,:)).', 'd', da);
hc = 3;                                                                % coarse far field of the same SDF (robot clearance
[xc, yc, zc] = meshgrid(single(-150:hc:210), single(-250:hc:150), single(-175:hc:240));   % checks beyond the fine grid)
dc = headSdf(xc, yc, zc, hc, mx);
% Collision distance = a guaranteed LOWER BOUND of the true distance to the skin.
% The blended SDF is accurate only close to the skin (blends and height fields
% distort it further away), so beyond 4 mm the exact Euclidean distance
% transform of the occupied voxels is used, minus half a voxel diagonal
% (voxel quantisation) and minus another half diagonal (trilinear interpolation).
hd = sqrt(3)*hc/2;
Lb = sqrt(edt3(dc <= 0, hc)) - 2*hd;
far = dc; far(dc >= 4) = max(Lb(dc >= 4), 0);
S.far = struct('x', squeeze(xc(1,:,1)), 'y', squeeze(yc(:,1,1)).', 'z', squeeze(zc(1,1,:)).', 'd', far);
Fe = griddedInterpolant({S.far.y, S.far.x, S.far.z}, double(Lb));
[xq, yq, zq] = meshgrid(S.x, S.y, S.z);
dfar = max(single(Fe(yq, xq, zq)), 0);
S.d(S.d >= 4) = dfar(S.d >= 4);
clear xq yq zq Fe Lb dc far dfar
cache = struct('F', F, 'V', V, 'N', N, 'C', C, 'B', B, 'S', S);
end

%% =============================== helpers ================================
function [F, V] = largestPart(F, V)
% keep only the largest connected component of a triangle mesh
n = size(V,1);
G = graph(F(:,[1 2 3]), F(:,[2 3 1]), [], n);
c = conncomp(G); keep = c == mode(c(F(:,1)));
fk = all(keep(F), 2); F = F(fk,:);
idx = zeros(n,1); idx(keep) = 1:nnz(keep);
V = V(keep,:); F = idx(F);
end

function V = taubin(F, V, w, iters)
% Taubin lambda|mu smoothing (no shrinkage), weighted per vertex by w in [0,1]
n = size(V,1);
A = sparse(F(:,[1 2 3 2 3 1]), F(:,[2 3 1 1 2 3]), 1, n, n); A = double(A > 0);
deg = max(full(sum(A, 2)), 1);
for k = 1:iters
    for f = [0.5 -0.53]
        V = V + f*w.*((A*V)./deg - V);
    end
end
end

function ao = ambientOcc(x, y, z, d, V, N, steps)
% cheap ambient occlusion: how far the SDF grows when stepping out along the normal
ao = zeros(size(V,1), 1);
for s = steps
    P = V + s*N;
    ds = double(interp3(x, y, z, d, single(P(:,1)), single(P(:,2)), single(P(:,3)), 'linear', s));
    ao = ao + min(1, max(0, ds/s));
end
ao = ao/numel(steps);
end

function d = vesica(x, y, a, b)
% 2-D SDF of an almond (intersection of two discs) with half-width a, half-height b
r = (a^2/b + b)/2; c = r - b;
d = max(sqrt(x.^2 + (y + c).^2) - r, sqrt(x.^2 + (y - c).^2) - r);
end

function d = sdEll(x,y,z,c,r)
k0 = sqrt(((x-c(1))/r(1)).^2 + ((y-c(2))/r(2)).^2 + ((z-c(3))/r(3)).^2);
d = (k0 - 1)*min(r);
end

function d = sdSup(x,y,z,c,r,p)
% superellipsoid (flatter than an ellipsoid for p > 2)
k0 = (abs((x-c(1))/r(1)).^p + abs((y-c(2))/r(2)).^p + abs((z-c(3))/r(3)).^p).^(1/p);
d = (k0 - 1)*min(r);
end

function d = sdCaps(x,y,z,a,b,r)
ba = b - a; px = x - a(1); py = y - a(2); pz = z - a(3);
t = min(max((px*ba(1) + py*ba(2) + pz*ba(3))/(ba*ba.'), 0), 1);
d = sqrt((px - t*ba(1)).^2 + (py - t*ba(2)).^2 + (pz - t*ba(3)).^2) - r;
end

function d = sdCone(x,y,z,a,b,r1,r2)
% capsule whose radius varies linearly from r1 (at a) to r2 (at b)
ba = b - a; px = x - a(1); py = y - a(2); pz = z - a(3);
t = min(max((px*ba(1) + py*ba(2) + pz*ba(3))/(ba*ba.'), 0), 1);
d = sqrt((px - t*ba(1)).^2 + (py - t*ba(2)).^2 + (pz - t*ba(3)).^2) - (r1 + t*(r2 - r1));
end

function d = smin(a, b, k)
hh = max(k - abs(a - b), 0)/k;
d = min(a, b) - hh.^2*k/4;
end

function d = smax(a, b, k)
d = -smin(-a, -b, k);
end

function hF = faceHeight(u, v)
% Height z (mm) of the front of the face over (u, v): u = lateral distance
% from the mid-line, v = yE (eye level 0, + towards the forehead).
au = abs(u);
pc = @(xn, yn, q) pchip(xn, yn, min(max(q, xn(1)), xn(end)));
G  = @(cu, cv, su, sv, p) exp(-(abs((u - cu)/su).^p + abs((v - cv)/sv).^p));
Gs = @(cu, cv, su, sv, p) G(cu, cv, su, sv, p) + G(-cu, cv, su, sv, p);
sig = @(t) 1./(1 + exp(-t));
% mid-line profile without nose and lips (forehead -> chin -> under the chin)
vP = [-128 -121 -113 -105 -96 -86 -74 -66 -58 -48 -40 -25 -10   0   6  18  30  45  60  75  90 105 120 136];
zP = [ -15    2   12   16  15  13  19  20  21  20  18  15  14  15 15.5 19.5 19 16 12   4 -12 -27 -52 -80];
Pm = pc(vP, zP, v);
% lateral fall-off: the face is rounder at the mouth and chin, flatter across the eyes
vA = [-124  -108   -90   -66   -50   -25     0     30     60    90    120];
aA = [0.020 0.012 0.010 0.012 0.009 0.005 0.003 0.0088 0.0088 0.0088 0.0088];
bA = [2e-5 9e-6 4e-6 1e-6 8.7e-7 1.56e-6 1.5e-6 6.8e-7 6.8e-7 6.8e-7 6.8e-7];
hF = Pm - pc(vA, aA, v).*u.^2 - pc(vA, bA, v).*u.^4;
% ---- forehead, brow ridge, temples ----
browL = exp(-((v - 15.5 - 0.0025*u.^2)/6.5).^2) .* exp(-(u/54).^6);
hF = hF + 3.0*browL.*(1 - 0.85*exp(-(u/14).^2)) + 1.2*G(0, 14, 9, 6, 2);   % two brow ridges + glabella
hF = hF - 3.0*Gs(61, 16, 9, 14, 2);
% ---- orbits (eyes sit in a hollow), tear trough, cheekbones, cheeks ----
hF = hF - 8*Gs(32, 1, 17, 12, 2.2);
hF = hF - 1.0*Gs(25, -12, 9, 2.5, 2);
hF = hF + 5.0*Gs(42, -20, 13, 11, 2);
hF = hF + 2.0*Gs(40, -52, 15, 15, 2);
% ---- nose: dorsum (variable width), tip lobule, alae with alar groove ----
vN = [-50 -47 -44 -40 -36 -30 -20 -10   0   6  12];
nN = [  0   6  14  19  20  19  15  10   5 2.5   0];
wN = pc([-50 -42 -36 -28 -15 0 16], [8 9.5 8.5 8 7 6.5 7.5], v);
nose = pc(vN, nN, v) .* exp(-0.69*(au./wN).^2.6);
nose = nose + 2.5*G(0, -38, 7, 5.5, 2);
for s = [-1 1]
    ra = sqrt(((u - s*13)/6.8).^2 + ((v + 41.5)/8.0).^2);
    ala = 5.2*exp(-ra.^2.2);
    nose = (nose.^6 + ala.^6).^(1/6);                              % soft max: ala merges with the tip
    nose = nose - 1.1*exp(-((ra - 1.05)/0.14).^2) .* sig((s*u - 13 + 4)/1.5);   % alar-facial groove
end
hF = hF + nose;
% ---- lips: philtrum, vermilion with cupid's bow, stomion, commissures ----
c  = max(1 - (u/25).^2, 0).^1.1;
c2 = max(1 - (u/23).^2, 0).^1.0;
vbU = -60 + 1.0*exp(-((au - 5)/3).^2) - 0.6*exp(-(u/1.8).^2) - 0.006*u.^2;   % upper vermilion border (cupid's bow)
hF = hF + 3.6*c .* exp(-((v + 63.4 + 0.003*u.^2)./(3.2 - 0.002*u.^2)).^2);   % upper lip
hF = hF + 5.0*c2 .* exp(-((v + 72.5 - 0.008*u.^2)./(4.2 - 0.004*u.^2)).^2);   % lower lip (fuller centre)
hF = hF + 0.6*c .* exp(-((v - vbU)/0.9).^2);                       % white roll
hF = hF - 3.0*c .* exp(-((v + 67)/0.9).^2);                        % stomion (mouth line)
hF = hF - 2.5*Gs(25, -67, 4, 3, 2);                                % commissures
hF = hF + 1.0*Gs(5, -54, 1.6, 5, 2) - 0.5*G(0, -54, 2, 5, 2);      % philtral columns + groove
% ---- nasolabial folds: cheek side slightly raised along a line ala -> past the commissure ----
for s = [-1 1]
    A = [s*18 -42]; Bp = [s*29 -73]; e = Bp - A; Le = norm(e); e = e/Le;
    t = ((u - A(1))*e(1) + (v - A(2))*e(2))/Le;
    sd = s*((u - A(1))*e(2) - (v - A(2))*e(1));                    % > 0 on the cheek side
    win = sig(t/0.08) .* sig((1 - t)/0.12);
    hF = hF + 1.3*sig(sd/1.8) .* win - 0.35*exp(-(sd/1.4).^2) .* win;
end
% ---- labiomental fold and chin ----
hF = hF - 2.5*exp(-((v + 87)/3.5).^2) .* max(1 - (u/20).^2, 0);
hF = hF + 3.0*G(0, -104, 15, 9, 2);
end

function [lu, ll] = lipMasks(u, v)
% soft masks of the upper / lower vermilion (same geometry as faceHeight)
au = abs(u); sig = @(t) 1./(1 + exp(-t));
vbU = -60 + 1.0*exp(-((au - 5)/3).^2) - 0.6*exp(-(u/1.8).^2) - 0.006*u.^2;
lu = sig((vbU - v)/0.45) .* sig((v + 67)/0.45) .* sig((24 - au)/0.8);
vbL = -77.5 + 0.011*u.^2;
ll = sig((v - vbL)/0.7) .* sig((-67 - v)/0.45) .* sig((23 - au)/0.8);
end

function d = headSdf(x, y, z, h, mx)
% signed distance (mm) of the head (face height field + volume + cranium,
% nostrils, ears, neck, eyelids, cap) on the grid x, y, z with spacing h
u = x - mx;
% ---- front surface of the face (height field) ----
U2 = double(u(:,:,1)); V2 = double(y(:,:,1));
hF = faceHeight(U2, V2);
[gx, gy] = gradient(hF, h);
dFront = (z - single(hF)) ./ single(sqrt(1 + gx.^2 + gy.^2));
% ---- face volume: outline seen from the front, jaw underside, rounded back ----
pc = @(xn, yn, q) pchip(xn, yn, min(max(q, xn(1)), xn(end)));
W2 = single(pc([-125 -117 -107 -95 -80 -60 -40 -20 0 20 50 80 100 115], ...
               [  10   30   43  53  59  64  68  70 71 70 66 62  52  35], V2));
zc = -30; Dz = 78;
k0 = ((abs(u)./W2).^3 + (abs(z - zc)/Dz).^3).^(1/3);
dBlk = (k0 - 1).*min(W2, Dz);
zv = squeeze(z(1,1,:));
vb = reshape(single(pc([-172 -100 -70 -50 -30 -10 10 40 60], [-90 -90 -96 -103 -110 -116 -120 -123 -123], double(zv))), 1, 1, []);
dBlk = smax(dBlk, vb - y, 7);                                     % lower border of the mandible
dBlk = smax(dBlk, y - 112, 25);
dFace = smax(dBlk, dFront, 9);
clear dBlk dFront k0
d = sdEll(x,y,z, [mx 28 -84], [73 98 88]);                        % cranium (behind the face)
d = smin(d, dFace, 10);
clear dFace
% ---- nostrils (on the under side of the nose) ----
for s = [-1 1]
    d = smax(d, -sdEll(x,y,z, [mx+s*5.4 -46.6 26.5], [3.3 2.0 3.3]), 1.0);
end
% ---- ears: plate + helix rim + antihelix + concha, back edge stands off the head ----
hel = [14 10; 24 4; 28 -6; 22 -15; 8 -18; -8 -16; -20 -11; -26 -5];     % helix rim (dy, dz) from the ear centre
ah  = [-14 -8; -4 -9; 8 -7; 17 -3; 20 4];                                % antihelix
hel = interp1(1:8, hel, linspace(1, 8, 22), 'spline');                     % smooth curves
ah  = interp1(1:5, ah, linspace(1, 5, 12), 'spline');
for s = [-1 1]
    ec = [mx+s*73 -12 -52];
    xe = x - s*0.28*max(ec(3) + 10 - z, 0);                                    % shear: back of the pinna flares out
    de = sdEll(xe,y,z, ec + [s*1 0 -3], [4.0 28 17]);                         % pinna plate
    for k = 1:size(hel,1)-1
        de = smin(de, sdCaps(xe,y,z, ec + [s*2.6 hel(k,:)], ec + [s*2.6 hel(k+1,:)], 3.0), 1.5);
    end
    de = smax(de, -sdEll(xe,y,z, ec + [s*6.2 2 -4], [3.0 19 10]), 1.5);     % scapha
    for k = 1:size(ah,1)-1
        de = smin(de, sdCaps(xe,y,z, ec + [s*2.4 ah(k,:)], ec + [s*2.4 ah(k+1,:)], 2.0), 1.5);
    end
    de = smax(de, -sdEll(xe,y,z, ec + [s*5.5 -6 3], [4.5 8.5 7]), 1.5);     % concha
    de = smin(de, sdEll(xe,y,z, ec + [s*0.5 -27 -4], [4.0 8 7]), 2);        % lobe
    de = smin(de, sdEll(xe,y,z, ec + [s*2.5 -6 11], [2.8 4 3]), 1.5);       % tragus
    de = smax(de, -sdEll(xe,y,z, ec + [s*3.5 -6 6], [3 2.5 2.5]), 1);       % ear canal
    de = smin(de, sdCaps(x,y,z, ec + [-s*9 -4 7], ec + [s*1 -2 6], 6), 3);   % ear root (crus) into the head
    d = smin(d, de, 3);
end
clear xe de
d = smin(d, sdCaps(x,y,z, [mx -98 -80], [mx -225 -102], 44), 24); % neck (ends under the blanket)
% ---- operated (right) eye: lid mass over the globe, fissure held open by the speculum ----
rr = sqrt(x.^2 + y.^2 + z.^2);
d = smin(d, rr - 15.0, 7);
ap = vesica(x, y, 13.0, 9.5);                                     % almond-shaped palpebral fissure
z0 = 9*(abs(x)/13).^4;                                            % shallower cut towards the canthi: the lids
d = smax(d, -max(ap, z0 - z), 1.8);                               % meet there, no gap beside the globe
d = smax(d, -(rr - 12.6), 1.0);                                   % room for the globe (R = 12 mm)
% ---- fellow (left) eye: closed lids, curved lid line, upper-lid crease ----
r2 = sqrt((x + 64).^2 + y.^2 + z.^2);
d = smin(d, r2 - 14.6, 7);
yl = -2 + 0.012*(x + 64).^2;
g = max(max(abs(y - yl) - 0.45, abs(x + 64) - 12.5), 7 - z);
d = smax(d, -g, 0.5);
yc = 7.5 - 0.010*(x + 64).^2;
g = max(max(abs(y - yc) - 0.5, abs(x + 64) - 12), 4 - z);
d = smax(d, -(g + 0.3), 0.8);
clear rr r2 ap g yl yc
% ---- surgical cap: 1.5 mm fabric over the cranium, soft wrinkles ----
cap = 1./(1 + exp(-(y - 66)/1.5));
d = d - cap.*(1.5 + 0.5*sin(x/7 + z/9).*sin(z/6 - x/11));
clear cap
d = smax(d, -(z + 169), 2);                                       % rests on the head rest
% normalise so |grad d| <= 1 (a distance that is never over-estimated -> safe clearance checks)
[ga, gb_, gc] = gradient(d, h);
d = d ./ max(sqrt(ga.^2 + gb_.^2 + gc.^2), 1);
end

function D = edt3(inside, h)
% exact squared Euclidean distance transform (mm^2) to the nearest occupied
% voxel: separable 1-D min-plus passes along each axis (brute force per line)
G = inf(size(inside)); G(inside) = 0;
for dim = 1:3
    perm = [dim setdiff(1:3, dim)];
    A = permute(G, perm); sz = size(A); n = sz(1); A = reshape(A, n, []);
    W = (h*((1:n).' - (1:n))).^2;                  % (i - j)^2 h^2
    out = zeros(size(A));
    for c0 = 1:400:size(A,2)
        c = c0:min(c0 + 399, size(A,2));
        out(:,c) = squeeze(min(reshape(A(:,c), 1, n, []) + W, [], 2));
    end
    G = ipermute(reshape(out, sz), perm);
end
D = G;
end
