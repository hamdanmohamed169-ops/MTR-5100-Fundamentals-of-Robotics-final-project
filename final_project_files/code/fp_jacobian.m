function [J, Jrcm] = fp_jacobian(q, P, dDH)
%FP_JACOBIAN  Geometric 6x3 Jacobian of the needle tip (Module 6).
%   J = [Jv; Jw] maps qdot to [v_tip; omega] in {0}.
%   Jrcm (3x3) maps qdot to the velocity of the needle material point that
%   is currently located at the nominal RCM (scleral entry, origin of {0}):
%       v_rcm = v_tip + omega x (0 - p) = (Jv + [p]x Jw) qdot.
%   For a perfect RCM mechanism this velocity is purely axial (q3dot*u).
n = size(P.DH, 1);
if nargin < 3, dDH = zeros(n,4); end
[T0n, Tall] = fp_fk(q, P, dDH);
p = T0n(1:3,4);
o = zeros(3,n+1); z = zeros(3,n+1);
z(:,1) = [0;0;1];
for i = 1:n
    o(:,i+1) = Tall(1:3,4,i);
    z(:,i+1) = Tall(1:3,3,i);
end
J = zeros(6,n);
for i = 1:n
    if P.DH(i,5) == 0          % revolute: [z x (p - o); z]
        J(1:3,i) = cross(z(:,i), p - o(:,i));
        J(4:6,i) = z(:,i);
    else                       % prismatic: [z; 0]
        J(1:3,i) = z(:,i);
    end
end
Sp   = [0 -p(3) p(2); p(3) 0 -p(1); -p(2) p(1) 0];
Jrcm = J(1:3,:) + Sp*J(4:6,:);
end
