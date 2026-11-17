function [M, c, G, V] = fp_arm_dynamics(qa, qad, P, withRotor)
%FP_ARM_DYNAMICS  Rigid-body dynamics of the 4-DOF positioning arm carrying
%   the RCM injection unit (Modules 2 and 6):
%       M(q) qdd + c(q,qd) + G(q) = tau
%   Bodies: arm links 1-4 (P.arm.link, COM and inertia in their DH frames)
%   and the RCM unit (P.arm.payload, rigidly fixed to frame {4}).
%   M = sum_b ( m_b Jv_b' Jv_b + Jw_b' R_b I_b R_b' Jw_b ) + diag(Jm N^2)
%   c_k = sum_ij ( dM_kj/dq_i - 1/2 dM_ij/dq_k ) qd_i qd_j   (Christoffel,
%         derivatives of M by central differences)
%   G = dV/dq,  V = -sum_b m_b g' p_cb,  g = [0 0 -9.81] in {E} (zE is up).
%   withRotor (default true) adds the reflected rotor inertia Jm N^2.
if nargin < 4, withRotor = true; end
[M, G, V] = massgrav(qa, P);
if withRotor, M = M + diag(P.arm.act.Jm .* P.arm.act.N.^2); end
if nargout > 1
    if nargin < 2 || isempty(qad), qad = zeros(4,1); end
    h = 1e-6; n = 4; dM = zeros(n,n,n);
    for i = 1:n
        e = zeros(n,1); e(i) = h;
        dM(:,:,i) = (massgrav(qa + e, P) - massgrav(qa - e, P))/(2*h);
    end
    Mdot = zeros(n); quad = zeros(n,1);
    for i = 1:n, Mdot = Mdot + dM(:,:,i)*qad(i); end
    for k = 1:n, quad(k) = qad.'*dM(:,:,k)*qad; end
    c = Mdot*qad - 0.5*quad;
end
end

function [M, G, V] = massgrav(qa, P)
[~, Ta] = fp_arm_fk(qa, P);
g = [0; 0; -P.g];
o = [P.arm.T_Eb(1:3,4) squeeze(Ta(1:3,4,1:3))];          % joint origins (axis i through o(:,i))
z = [P.arm.T_Eb(1:3,3) squeeze(Ta(1:3,3,1:3))];          % joint axes
bodies = [P.arm.link, P.arm.payload];
frameOf = [1 2 3 4 4];
M = zeros(4); G = zeros(4,1); V = 0;
for b = 1:numel(bodies)
    f = frameOf(b); R = Ta(1:3,1:3,f);
    pc = Ta(1:3,4,f) + R*bodies(b).rc;
    Jv = zeros(3,4); Jw = zeros(3,4);
    for j = 1:f
        Jv(:,j) = cross(z(:,j), pc - o(:,j)); Jw(:,j) = z(:,j);
    end
    m = bodies(b).m; I = bodies(b).I;
    M = M + m*(Jv.'*Jv) + Jw.'*(R*I*R.')*Jw;
    G = G - m*Jv.'*g;
    V = V - m*g.'*pc;
end
M = 0.5*(M + M.');
end
