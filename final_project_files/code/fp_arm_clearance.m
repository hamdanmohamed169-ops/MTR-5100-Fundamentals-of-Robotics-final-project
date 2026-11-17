function [dmin, info] = fp_arm_clearance(qa, P, CM, ENV, full, q)
%FP_ARM_CLEARANCE  Smallest distance (mm) between the positioning arm + the
%   RCM unit it carries (RCM robot at q_home) and the obstacles of FP_ARM_ENV,
%   for arm joints qa (4x1, rad); optional q = RCM joints (default q_home).
%   full = true uses every mesh vertex (verification), otherwise the voxel-
%   sampled planning model of FP_ARM_COLLISION_MODEL.
%   info.body : 1-4 arm links, 5 RCM unit;  info.obst : obstacle id;
%   info.p    : closest robot point in {E} (mm)
if nargin < 5 || isempty(full), full = false; end
[TE4, Ta] = fp_arm_fk(qa, P);
if full, L = CM.linkFull; U = CM.unitFull; else, L = CM.link; U = CM.unit; end
pts = cell(1,5);
for i = 1:4
    pts{i} = Ta(1:3,1:3,i)*L{i} + Ta(1:3,4,i)*1e3;
end
if nargin < 6 || isempty(q)
    pts{5} = TE4(1:3,1:3)*U + TE4(1:3,4)*1e3;
else                                                % RCM robot away from home: move its bodies
    if full, B = CM.ubFull; else, B = CM.ub; end
    TE0 = TE4*P.arm.T_40; [~, Tq] = fp_fk(q, P);
    Tb = {TE0, TE0*Tq(:,:,1), TE0*Tq(:,:,2), TE0*Tq(:,:,3)}; pu = cell(1,4);
    for i = 1:4, pu{i} = Tb{i}(1:3,1:3)*B{i} + Tb{i}(1:3,4)*1e3; end
    pts{5} = [pu{:}];
end
% needle points inside the eyeball are not checked here: inside the globe the
% in-eye safety layer applies (lens / retina / keep-out clearances, FP_LENS_CHECK)
rG = P.eye.R*1e3 + 0.6;
for i = 1:5, pts{i} = pts{i}(:, sum(pts{i}.^2, 1) > rG^2); end
n = cellfun(@(c) size(c,2), pts);
Pall = [pts{:}];
[d, id] = ENV.sdf(Pall);
[dmin, k] = min(d);
body = repelem(1:5, n);
info = struct('body', body(k), 'obst', id(k), 'p', Pall(:,k));
end
