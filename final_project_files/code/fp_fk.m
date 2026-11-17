function [T0n, Tall] = fp_fk(q, P, dDH)
%FP_FK  DH forward kinematics of the 3-DOF RCM robot (Module 4/5).
%   [T03, Tall] = fp_fk(q, P) returns the tool (needle-tip) pose T03 in {0}
%   and Tall(:,:,i) = T0i for i = 1..3.
%   dDH (optional, 3x4) adds errors to the DH rows [theta d a alpha]; it is
%   used to study manufacturing / assembly errors and RCM drift.
n = size(P.DH, 1);
if nargin < 3, dDH = zeros(n,4); end
DH = P.DH(:,1:4) + dDH;
type = P.DH(:,5);
T = eye(4);
Tall = zeros(4,4,n);
for i = 1:n
    th = DH(i,1); d = DH(i,2);
    if type(i) == 0, th = th + q(i); else, d = d + q(i); end
    T = T * dh_std(th, d, DH(i,3), DH(i,4));
    Tall(:,:,i) = T;
end
T0n = T;
end

function A = dh_std(th, d, a, al)
% Standard DH link transform A = Rz(th) Tz(d) Tx(a) Rx(al)
ct = cos(th); st = sin(th); ca = cos(al); sa = sin(al);
A = [ct -st*ca  st*sa a*ct;
     st  ct*ca -ct*sa a*st;
     0   sa     ca    d;
     0   0      0     1];
end
