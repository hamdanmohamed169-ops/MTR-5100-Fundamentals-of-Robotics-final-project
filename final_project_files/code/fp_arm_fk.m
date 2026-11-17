function [T, Tall] = fp_arm_fk(qa, P)
%FP_ARM_FK  Forward kinematics of the 4-DOF positioning arm (standard DH).
%   T = flange pose {4} in the eye frame {E} (metres); Tall(:,:,i) = pose
%   of arm frame i in {E} (i = 1..4). P.arm.DH rows = [theta_off d a alpha].
DH = P.arm.DH; n = size(DH,1); T = P.arm.T_Eb; Tall = zeros(4,4,n);
for i = 1:n
    th = DH(i,1) + qa(i); d = DH(i,2); a = DH(i,3); al = DH(i,4);
    ct = cos(th); st = sin(th); ca = cos(al); sa = sin(al);
    T = T*[ct -st*ca st*sa a*ct; st ct*ca -ct*sa a*st; 0 sa ca d; 0 0 0 1];
    Tall(:,:,i) = T;
end
end
