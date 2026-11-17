function TE0 = fp_scene_update(H, q, P, thetaEye, Vinj, qa)
%FP_SCENE_UPDATE  Move the prototype scene to RCM joints q = [q1; q2; q3]
%   (rad, rad, m) and, optionally, positioning-arm joints qa (4x1, rad).
%   thetaEye : eye rotation about xE (rad) for the twitch event
%   Vinj     : injected volume (uL) -> plunger displacement
%   qa       : arm joints; the RCM base {0} then rides on the arm flange,
%              T_E0 = T_E4(qa) * T_40. Returns that T_E0 (metres).
if nargin < 4, thetaEye = 0; end
if nargin < 5, Vinj = 0; end
mmT = @(T) [T(1:3,1:3) T(1:3,4)*1e3; 0 0 0 1];
TE0 = P.T_E0;
if nargin >= 6 && ~isempty(qa)
    [TE4, Ta] = fp_arm_fk(qa, P);
    TE0 = TE4*P.arm.T_40;
    if isfield(H, 'arm') && ~isempty(H.arm)
        for i = 1:4, H.arm.T(i+1).Matrix = mmT(Ta(:,:,i)); end
    end
end
H.base.Matrix = mmT(TE0);
[~, Tall] = fp_fk(q, P);
H.L1.Matrix = mmT(TE0*Tall(:,:,1)); H.L2.Matrix = mmT(TE0*Tall(:,:,2));
H.L3.Matrix = mmT(TE0*Tall(:,:,3));
if isfield(H, 'cable') && isgraphics(H.cable)               % flexible cable: tool unit -> J1 housing port
    T2 = mmT(TE0*Tall(:,:,2)); TB = mmT(TE0);
    p0 = T2*[0; 0; -128; 1]; d0 = -T2(1:3,3);
    p3 = TB*[H.cablePort; 1]; d3 = TB(1:3,1:3)*H.cableDir;
    pts = H.bez(p0(1:3), p0(1:3) + 45*d0, p3(1:3) + 45*d3, p3(1:3), 32);
    [X, Y, Z] = H.tubeMesh(pts, 2.3);
    set(H.cable, 'XData', X, 'YData', Y, 'ZData', Z);
end
A = pi*(P.iop.syr_D*1e3)^2/4;                               % mm^2
H.plunger.Matrix = makehgtform('translate', [0 0 Vinj/A]);  % uL = mm^3
c = cos(thetaEye); s = sin(thetaEye);
H.eye.Matrix = [1 0 0 0; 0 c -s 0; 0 s c 0; 0 0 0 1];
end
