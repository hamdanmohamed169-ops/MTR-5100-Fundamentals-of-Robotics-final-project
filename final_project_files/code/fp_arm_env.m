function ENV = fp_arm_env(P)
%FP_ARM_ENV  Obstacles around the positioning arm, as signed distance (mm, {E}):
%     1 patient head + neck   (fine + coarse SDF grids of FP_HEAD_MESH, gradient-normalised)
%     2 patient body + blanket (SDF grid of FP_HEAD_MESH)
%     3 operating table, head rest and pedestal (boxes)
%     4 stereo-camera stand: column, boom, camera head (capsules + box)
%   ENV.sdf(p) -> [d, id] for points p (3xN, mm): distance to the nearest
%   obstacle and which one (1..4). Cached per MATLAB session.
persistent cache
if ~isempty(cache), ENV = cache; return; end
[~, Vh, ~, ~, Bm, S] = fp_head_mesh();
Fh = griddedInterpolant({S.y, S.x, S.z}, double(S.d), 'linear', 'none');
Ff = griddedInterpolant({S.far.y, S.far.x, S.far.z}, double(S.far.d), 'linear', 'none');
Fb = griddedInterpolant({S.body.y, S.body.x, S.body.z}, double(S.body.d), 'linear', 'none');
gH = [S.x([1 end]); S.y([1 end]); S.z([1 end])];
gF = [S.far.x([1 end]); S.far.y([1 end]); S.far.z([1 end])];
gB = [S.body.x([1 end]); S.body.y([1 end]); S.body.z([1 end])];
bbH = [min(Vh, [], 1).' max(Vh, [], 1).'];          % bounding boxes of the meshes (far-away lower bound)
bbB = [min(Bm.V, [], 1).' max(Bm.V, [], 1).'];
% table: mattress + frame + head rest + pedestal (axis-aligned boxes [lo hi])
boxes = {[-345 -440 -198; 175 210 -170], [-335 -432 -220; 165 202 -196], ...
         [-66 -70 -218; 2 70 -150], [-150 -260 -560; -20 -120 -196]};
% camera stand (same geometry as FP_SCENE / buildCamera)
mm = 1e3; K = P.cam; Cc = K.C*mm; tg = K.target*mm;
zc = (tg - Cc)/norm(tg - Cc); bx = cross(zc, [0;0;1]); bx = bx/norm(bx); by = cross(zc, bx);
pole = [-300; 300]; zf = -560;
M = Cc - 125*zc; J0 = [pole; 640]; J1 = [M(1:2) + 0.45*(pole - M(1:2)); 700];
caps = {[pole; zf+40], [pole; 640], 34; J0, J1 + [0;0;10], 30; J1, M + [0;0;30], 30; M + [0;0;30], Cc, 22};
camR = [bx by zc]; camLo = [-98; -34; -75]; camHi = [98; 34; 23];
ENV.sdf = @envSdf;
ENV.names = {'patient head', 'patient body', 'table', 'camera stand'};
cache = ENV;

    function [d, id] = envSdf(p)
        n = size(p, 2); D = inf(4, n);
        D(1,:) = gridSdf({Fh, Ff}, {gH, gF}, bbH, p);
        D(2,:) = gridSdf({Fb}, {gB}, bbB, p);
        for b = 1:numel(boxes), D(3,:) = min(D(3,:), boxSdf(p, boxes{b}(1,:).', boxes{b}(2,:).')); end
        for c = 1:size(caps,1), D(4,:) = min(D(4,:), capSdf(p, caps{c,1}, caps{c,2}, caps{c,3})); end
        D(4,:) = min(D(4,:), boxSdf(camR.'*(p - Cc), camLo, camHi));
        [d, id] = min(D, [], 1);
    end
end

function d = gridSdf(Fs, gs, bb, p)
% SDF from nested grids (fine first); outside every grid the distance to the
% mesh bounding box bb is used (a lower bound of the true distance)
d = nan(1, size(p,2)); todo = true(1, size(p,2));
for i = 1:numel(Fs)
    g = gs{i}; in = todo & all(p >= g(:,1) & p <= g(:,2), 1);
    d(in) = Fs{i}(p(2,in), p(1,in), p(3,in)); todo(in) = false;
end
if any(todo)
    q = p(:,todo); e = max(max(bb(:,1) - q, q - bb(:,2)), 0);
    d(todo) = sqrt(sum(e.^2, 1));
end
end

function d = boxSdf(p, lo, hi)
c = (lo + hi)/2; h = (hi - lo)/2;
q = abs(p - c) - h;
d = sqrt(sum(max(q, 0).^2, 1)) + min(max(q, [], 1), 0);
end

function d = capSdf(p, a, b, r)
ba = b - a; t = min(max((ba.'*(p - a))/(ba.'*ba), 0), 1);
d = sqrt(sum((p - a - ba*t).^2, 1)) - r;
end
