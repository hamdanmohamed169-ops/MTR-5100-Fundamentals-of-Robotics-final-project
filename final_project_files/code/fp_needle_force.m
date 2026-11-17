function [F, st] = fp_needle_force(depth, vel, st, P)
%FP_NEEDLE_FORCE  Axial needle-tissue force model (resists insertion).
%   depth : tip position beyond the scleral surface (= q3, m)
%   vel   : insertion speed (m/s)
%   st    : state struct with field .punctured (logical)
%   Before puncture the sclera tents like a spring (F = k*depth). When F
%   reaches the puncture force the wall is pierced, the force drops and only
%   shaft friction + vitreous drag remain (classic force-drop signature that
%   the force sensor uses to detect entry into the vitreous).
N = P.needle;
if depth <= 0
    F = 0; st.punctured = false;   % needle outside the eye
    return
end
if ~st.punctured
    F = N.k_scl*depth;
    if F >= N.F_punct
        st.punctured = true;
        st.depth_punct = depth;
    end
end
if st.punctured
    F = N.f_fric*depth + N.c_vit*vel;
end
end
