function CM = fp_arm_collision_model(P)
%FP_ARM_COLLISION_MODEL  Point-cloud model of the moving bodies of the
%   4-DOF positioning arm and of the RCM injection unit it carries, taken
%   from the 3-D prototype geometry itself (FP_SCENE), so the planner checks
%   exactly what is drawn.
%   CM.link{i}  : 3xNi surface points of arm body i (i = 1..4) in its DH frame {i} (mm)
%   CM.unit     : 3xNu points of the RCM unit (robot at q_home) in the flange frame {4} (mm)
%   CM.ub{i}    : the RCM robot bodies (base, L1, L2, L3) in their frames {0}..{3} (any q)
%   CM.*Full    : the same, all mesh vertices (final verification)
%   CM.voxel    : sampling voxel size (mm) of the planning model, CM.sampErr its
%                 worst-case sampling error (half the voxel diagonal)
%   Cached per MATLAB session.
persistent cache
if ~isempty(cache), CM = cache; return; end
vox = [6 6 6 6]; voxU = 4;
fig = figure('Visible', 'off'); ax = axes(fig);
H = fp_scene(ax, P, struct('head', false, 'setup', true));
qa = deg2rad([20; 60; -80; -40]);                    % any reference configuration
mmT = @(T) [T(1:3,1:3) T(1:3,4)*1e3; 0 0 0 1];
fp_scene_update(H, P.q_home, P, 0, 0, qa);
[TE4, Ta] = fp_arm_fk(qa, P);
CM.link = cell(1,4); CM.linkFull = cell(1,4);
for i = 1:4
    W = collect(H.arm.T(i+1), eye(4));               % world points (mm)
    Pl = mmT(Ta(:,:,i)) \ [W; ones(1, size(W,2))];
    CM.linkFull{i} = Pl(1:3,:);
    CM.link{i} = voxelSample(Pl(1:3,:), vox(i));
end
hb = {H.base, H.L1, H.L2, H.L3};
TE0 = TE4*P.arm.T_40; [~, Tq] = fp_fk(P.q_home, P);
Tb = {TE0, TE0*Tq(:,:,1), TE0*Tq(:,:,2), TE0*Tq(:,:,3)};
CM.ub = cell(1,4); CM.ubFull = cell(1,4); W = zeros(3,0);
for i = 1:4                                          % RCM-robot bodies in their own frames {0},{1},{2},{3}
    Wi = collect(hb{i}, eye(4)); W = [W Wi]; %#ok<AGROW>
    Pi = mmT(Tb{i}) \ [Wi; ones(1, size(Wi,2))];
    CM.ubFull{i} = Pi(1:3,:); CM.ub{i} = voxelSample(Pi(1:3,:), voxU);
end
Pu = mmT(TE4) \ [W; ones(1, size(W,2))];
CM.unitFull = Pu(1:3,:);
CM.unit = voxelSample(Pu(1:3,:), voxU);
CM.voxel = [vox voxU]; CM.sampErr = sqrt(3)*max([vox voxU])/2;
CM.names = {'arm link 1', 'arm link 2', 'arm link 3', 'arm link 4', 'RCM unit'};
close(fig);
cache = CM;
end

function W = collect(h, M)
% all patch / surface vertices below graphics object h, in world coordinates
W = zeros(3,0);
if isa(h, 'matlab.graphics.primitive.Transform'), M = M*h.Matrix; end
switch class(h)
    case 'matlab.graphics.primitive.Patch'
        v = h.Vertices.'; v = v(:, all(isfinite(v), 1));
        W = M(1:3,1:3)*v + M(1:3,4);
    case 'matlab.graphics.chart.primitive.Surface'
        v = [h.XData(:) h.YData(:) h.ZData(:)].'; v = v(:, all(isfinite(v), 1));
        W = M(1:3,1:3)*v + M(1:3,4);
end
for c = h.Children.'
    W = [W collect(c, M)]; %#ok<AGROW>
end
end

function Q = voxelSample(Pts, s)
% one point (the mean) per occupied voxel of size s
[~, ~, g] = unique(round(Pts.'/s), 'rows');
Q = [accumarray(g, Pts(1,:).', [], @mean) accumarray(g, Pts(2,:).', [], @mean) accumarray(g, Pts(3,:).', [], @mean)].';
end
