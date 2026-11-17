function [M, c, G, V] = fp_dynamics(q, qd, P, scale)
%FP_DYNAMICS  Rigid-body dynamics  M(q) qdd + c(q,qd) + G(q) = tau  (Module 6).
%   M : 3x3 inertia matrix from COM Jacobians (Lagrange formulation)
%       M = sum_i ( m_i Jv_i' Jv_i + Jw_i' R_i I_i R_i' Jw_i )
%   c : Coriolis/centrifugal vector from Christoffel symbols, computed with
%       central differences of M:  c_k = sum_ij (dM_kj/dq_i - 1/2 dM_ij/dq_k) qd_i qd_j
%   G : gravity vector  G = dV/dq,  V = -sum_i m_i g0' p_ci
%   scale (optional) multiplies all masses/inertias -> model-uncertainty study.
if nargin < 4, scale = 1; end
[M, G, V] = massgrav(q, P, scale);
% Reflected rotor inertia (Module 2): J_m * N^2 seen at the joint side.
% Constant -> adds to M but not to the Coriolis terms.
M = M + diag(P.act.Jm .* P.act.Nr.^2);
if nargout > 1
    h = 1e-6; n = size(P.DH,1);
    dM = zeros(n,n,n);
    for i = 1:n
        e = zeros(n,1); e(i) = h;
        dM(:,:,i) = (massgrav(q+e,P,scale) - massgrav(q-e,P,scale))/(2*h);
    end
    Mdot = zeros(n);
    for i = 1:n, Mdot = Mdot + dM(:,:,i)*qd(i); end
    quad = zeros(n,1);
    for k = 1:n, quad(k) = qd.'*dM(:,:,k)*qd; end
    c = Mdot*qd - 0.5*quad;
end
end

function [M, G, V] = massgrav(q, P, scale)
[~, Tall] = fp_fk(q, P);
g0 = P.T_E0(1:3,1:3).' * P.eye.gE;        % gravity expressed in {0}
n = size(P.DH,1);
o = zeros(3,n+1); z = zeros(3,n+1); z(:,1) = [0;0;1];
for i = 1:n
    o(:,i+1) = Tall(1:3,4,i);  z(:,i+1) = Tall(1:3,3,i);
end
M = zeros(n); G = zeros(n,1); V = 0;
for b = 1:n
    R  = Tall(1:3,1:3,b);
    pc = Tall(1:3,4,b) + R*P.link(b).rc;
    Jv = zeros(3,n); Jw = zeros(3,n);
    for j = 1:b
        if P.DH(j,5) == 0
            Jv(:,j) = cross(z(:,j), pc - o(:,j));  Jw(:,j) = z(:,j);
        else
            Jv(:,j) = z(:,j);
        end
    end
    m = scale*P.link(b).m;  I = scale*P.link(b).I;
    M = M + m*(Jv.'*Jv) + Jw.'*(R*I*R.')*Jw;
    G = G - m*Jv.'*g0;
    V = V - m*g0.'*pc;
end
M = 0.5*(M + M.');
end
