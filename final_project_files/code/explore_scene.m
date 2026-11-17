function explore_scene()
%EXPLORE_SCENE  Walk around the 3-D operating room without running the procedure.
%   One click opens the room (patient, table, stereo camera, 4-DOF arm and the
%   3-DOF RCM injection unit) with:
%     - preset views (room, patient + robot, eye close-up, RCM unit, arm,
%       camera stand, surgeon's view, what the stereo camera sees)
%     - free navigation:  mouse  left-drag = orbit, right-drag = pan, wheel = zoom
%                         keys   W/S forward/back, A/D sideways, Q/E down/up,
%                                arrows = look around, R = reset view
%     - robot poses (parked, pre-dock, docked, needle at site A / B) and joint
%       sliders for both robots with a live clearance read-out to the patient
%       (arm + unit vs face, body, table, camera stand)
%     - "Plan path to the eye": RRT-Connect from the CURRENT arm pose to the
%       pre-dock pose, animated with its speed limit (FP_ARM_PLAN)
%     - "Play live demo" / "Record video" (DEMO_LIVE) and "Run full analysis"
%       (MAIN_FINAL_PROJECT)
here = fileparts(mfilename('fullpath')); if isempty(here), here = pwd; end
addpath(here);
P = fp_params(); mm = 1e3;
bg = [0.05 0.07 0.10]; pan = [0.09 0.11 0.15]; fg = [0.88 0.91 0.95]; acc = [0.35 0.80 1.0]; gold = [1 0.85 0.3];
fig = figure('Name', 'MTR 5100 - explore the operating room', 'Color', bg, 'Position', [30 40 1520 860], ...
    'NumberTitle', 'off', 'MenuBar', 'none', 'ToolBar', 'figure', 'InvertHardcopy', 'off');
wb = waitbar(0.1, 'Building the operating room (patient, robot, camera)...');
% ---------------- 3-D room ----------------
pM = uipanel(fig, 'Position', [0 0 0.76 1], 'BorderType', 'none', 'BackgroundColor', bg);
axBG = axes(pM, 'Position', [0 0 1 1]);
gy = linspace(0, 1, 200).'; gx = linspace(-1, 1, 300);
G = 0.04 + 0.10*(1 - gy).^1.5 .* (1 - 0.6*gx.^2);
image(axBG, cat(3, 0.55*G, 0.75*G, 1.0*G)); axis(axBG, 'off');
ax = axes(pM, 'Position', [0 0 1 1], 'Color', 'none');
H = fp_scene(ax, P);
axis(ax, 'off'); daspect(ax, [1 1 1]); ax.Clipping = 'off';
xlim(ax, [-1300 1100]); ylim(ax, [-1300 1100]); zlim(ax, [-570 900]);
light(ax, 'Position', [0.4 -1 2], 'Style', 'infinite', 'Color', [1 0.97 0.92]);
light(ax, 'Position', [-1 0.8 0.5], 'Style', 'infinite', 'Color', [0.45 0.5 0.6]);
light(ax, 'Position', [0.2 1 -0.3], 'Style', 'infinite', 'Color', [0.25 0.35 0.5]);
lighting(ax, 'gouraud'); ax.AmbientLightColor = [0.85 0.88 0.95];
waitbar(0.5, wb, 'Collision model and obstacle distance field...');
CM = fp_arm_collision_model(P); ENV = fp_arm_env(P);
hNear = plot3(ax, nan, nan, nan, 'o', 'MarkerSize', 9, 'LineWidth', 2, 'Color', [1 0.3 0.25]);
hPath = plot3(ax, nan, nan, nan, ':', 'Color', gold, 'LineWidth', 2);
hRaw  = plot3(ax, nan, nan, nan, '-', 'Color', [1 0.35 0.3 0.55], 'LineWidth', 1.2);
hTrail = animatedline(ax, 'Color', acc, 'LineWidth', 2.4);
% ---------------- key robot poses ----------------
TFdock = P.T_E0/P.arm.T_40;
qDock = fp_arm_ik(TFdock, P, zeros(4,1));
TFpre = TFdock; TFpre(1:3,4) = TFpre(1:3,4) - 0.040*P.T_E0(1:3,1);
qPre = fp_arm_ik(TFpre, P, qDock);
qPark = fp_arm_ik([0.34; -0.36; 0.05; -pi/2], P, deg2rad([60; 80; -100; -70]));
qA = [deg2rad(-20); deg2rad(100); 7e-3]; qB = [deg2rad(20); deg2rad(100); 7e-3];
poses = struct('name', {'Parked (beside the chest)', 'Pre-dock (40 mm before the eye)', 'Docked (RCM on the entry point)', ...
                        'Needle inserted - site A', 'Needle inserted - site B'}, ...
               'qa', {qPark, qPre, qDock, qDock, qDock}, 'q', {P.q_home, P.q_home, P.q_home, qA, qB});
st.qa = qPark; st.q = P.q_home; st.busy = false;
% ---------------- preset views {target, az, el, distance, view angle} ----------------
pE = P.eye.pE.'*mm;
views = {'Operating room (overview)', {[-100 -150 50], 50, 24, 4200, 34};
         'Patient + robot',           {[60 -80 150], 60, 22, 2800, 24};
         'Face and eye',              {[-30 -30 20], 30, 50, 1400, 16};
         'Eye close-up (speculum)',   {[0 0 8], 25, 55, 1400, 3.2};
         'RCM injection unit',        {[-10 -5 60], 140, 25, 1400, 12};
         'Needle entry (pars plana)', {pE + [-2 0 2], 150, 35, 1400, 4.5};
         '4-DOF arm (side)',          {[150 -60 150], 95, 8, 2800, 22};
         'Stereo camera stand',       {[-150 150 250], 210, 20, 2800, 26};
         "Surgeon's view (head end)", {[-20 0 20], 180, 62, 1400, 14};
         'Top view',                  {[0 -100 0], 0, 89, 3000, 30}};
% ---------------- control panel ----------------
pC = uipanel(fig, 'Position', [0.76 0 0.24 1], 'BackgroundColor', pan, 'BorderType', 'none');
y = 0.965;
lbl = @(txt, yy) uicontrol(pC, 'Style', 'text', 'Units', 'normalized', 'Position', [0.05 yy 0.9 0.025], ...
    'String', txt, 'BackgroundColor', pan, 'ForegroundColor', acc, 'FontWeight', 'bold', 'FontSize', 10, 'HorizontalAlignment', 'left');
lbl('VIEWS', y); y = y - 0.035;
for i = 1:size(views,1)
    c = mod(i-1, 2); r = floor((i-1)/2);
    uicontrol(pC, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.05 + 0.46*c, y - 0.036*r, 0.44, 0.032], ...
        'String', views{i,1}, 'FontSize', 8.5, 'Callback', @(~,~) goView(views{i,2}));
end
y = y - 0.036*ceil(size(views,1)/2) - 0.01;
uicontrol(pC, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.05 y 0.44 0.032], 'String', 'View from the stereo camera', ...
    'FontSize', 8.5, 'Callback', @(~,~) cameraView());
uicontrol(pC, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.51 y 0.44 0.032], 'String', 'Walk in from the door', ...
    'FontSize', 8.5, 'Callback', @(~,~) goView({[-100 -150 80], 35, 4, 5200, 50}));
y = y - 0.05;
lbl('ROBOT POSE', y); y = y - 0.035;
uicontrol(pC, 'Style', 'popupmenu', 'Units', 'normalized', 'Position', [0.05 y 0.9 0.03], 'String', {poses.name}, ...
    'FontSize', 9, 'Callback', @(s,~) setPose(s.Value));
y = y - 0.045;
lbl('JOINTS  (4-DOF arm | 3-DOF RCM robot)', y); y = y - 0.03;
jn = {'arm q1 base (deg)', 'arm q2 shoulder (deg)', 'arm q3 elbow (deg)', 'arm q4 wrist (deg)', ...
      'RCM q1 yaw (deg)', 'RCM q2 pitch (deg)', 'RCM q3 depth (mm)'};
lo = [rad2deg(P.arm.qmin); rad2deg(P.qmin(1:2)); P.qmin(3)*mm];
hi = [rad2deg(P.arm.qmax); rad2deg(P.qmax(1:2)); P.qmax(3)*mm];
sl = gobjects(1,7); sv = gobjects(1,7);
for j = 1:7
    uicontrol(pC, 'Style', 'text', 'Units', 'normalized', 'Position', [0.05 y 0.5 0.022], 'String', jn{j}, ...
        'BackgroundColor', pan, 'ForegroundColor', fg, 'FontSize', 8.5, 'HorizontalAlignment', 'left');
    sv(j) = uicontrol(pC, 'Style', 'text', 'Units', 'normalized', 'Position', [0.6 y 0.35 0.022], 'String', '', ...
        'BackgroundColor', pan, 'ForegroundColor', gold, 'FontSize', 8.5, 'HorizontalAlignment', 'right');
    sl(j) = uicontrol(pC, 'Style', 'slider', 'Units', 'normalized', 'Position', [0.05 y - 0.022 0.9 0.02], ...
        'Min', lo(j), 'Max', hi(j), 'Value', lo(j), 'Callback', @(~,~) fromSliders());
    addlistener(sl(j), 'ContinuousValueChange', @(~,~) fromSliders());
    y = y - 0.045;
end
txtClr = uicontrol(pC, 'Style', 'text', 'Units', 'normalized', 'Position', [0.05 y - 0.085 0.9 0.11], 'String', '', ...
    'BackgroundColor', [0.05 0.07 0.10], 'ForegroundColor', [0.4 1 0.5], 'FontName', 'Consolas', 'FontSize', 8.5, ...
    'HorizontalAlignment', 'left');
y = y - 0.12;
lbl('RUN', y); y = y - 0.04;
uicontrol(pC, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.05 y 0.9 0.036], ...
    'String', 'Plan path to the eye from this arm pose (RRT-Connect)', 'FontSize', 9, 'FontWeight', 'bold', ...
    'Callback', @(~,~) planHere());
y = y - 0.042;
uicontrol(pC, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.05 y 0.44 0.036], 'String', 'Play live demo', ...
    'FontSize', 9, 'FontWeight', 'bold', 'BackgroundColor', [0.25 0.6 0.35], 'ForegroundColor', 'w', 'Callback', @(~,~) demo_live());
uicontrol(pC, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.51 y 0.44 0.036], 'String', 'Record demo video', ...
    'FontSize', 9, 'Callback', @(~,~) demo_live('record'));
y = y - 0.042;
uicontrol(pC, 'Style', 'pushbutton', 'Units', 'normalized', 'Position', [0.05 y 0.9 0.036], ...
    'String', 'Run full engineering analysis (main_final_project)', 'FontSize', 9, 'Callback', @(~,~) runMain());
axO = axes(pM, 'Position', [0 0 1 1], 'Visible', 'off', 'XLim', [0 1], 'YLim', [0 1], 'HitTest', 'off', ...
    'PickableParts', 'none');                                    % 2-D overlay for the title and help line
text(axO, 0.008, 0.015, ['Mouse: left-drag orbit | right-drag pan | wheel zoom     ' ...
    'Keys: W/S forward-back  A/D sideways  Q/E down-up  arrows look around  R reset'], 'Color', [0.6 0.66 0.75], ...
    'FontSize', 9, 'VerticalAlignment', 'bottom');
ttl = text(axO, 0.012, 0.975, 'MTR 5100  |  Operating room - free exploration', 'Color', acc, 'FontSize', 14, ...
    'FontWeight', 'bold', 'VerticalAlignment', 'top', 'Interpreter', 'none');
% ---------------- navigation ----------------
nav.mode = ''; nav.p0 = [0 0];
set(fig, 'WindowButtonDownFcn', @mouseDown, 'WindowButtonUpFcn', @mouseUp, ...
    'WindowButtonMotionFcn', @mouseMove, 'WindowScrollWheelFcn', @wheel, 'KeyPressFcn', @key);
close(wb);
setPose(1);
goView(views{2,2}, 0);

%% ======================= callbacks ======================================
    function applyPose(extra)
        fp_scene_update(H, st.q, P, 0, 0, st.qa);
        v = [rad2deg(st.qa); rad2deg(st.q(1:2)); st.q(3)*mm];
        for jj = 1:7
            sl(jj).Value = min(max(v(jj), lo(jj)), hi(jj));
            sv(jj).String = sprintf('%.1f', v(jj));
        end
        [d, inf_] = fp_arm_clearance(st.qa, P, CM, ENV, false, st.q);
        [dS, infS] = fp_arm_selfclear(st.qa, P, CM);
        if dS < 10, d = min(d, dS); end
        if d >= 25, col = [0.4 1 0.5]; elseif d >= 10, col = [1 0.85 0.3]; else, col = [1 0.35 0.3]; end
        txtClr.ForegroundColor = col;
        if nargin < 1, extra = ''; end
        txtClr.String = sprintf(' clearance %6.1f mm  (lower bound)\n closest: %s -> %s\n self-collision %6.1f mm\n (%s)\n needle tip %s\n%s', ...
            d, CM.names{inf_.body}, ENV.names{inf_.obst}, dS, infS.pair, tipText(), extra);
        set(hNear, 'XData', inf_.p(1), 'YData', inf_.p(2), 'ZData', inf_.p(3));
        drawnow limitrate;
    end

    function s = tipText()
        if st.q(3) > 0, s = sprintf('%.1f mm inside the eye', st.q(3)*mm); else, s = sprintf('%.1f mm outside', -st.q(3)*mm); end
    end

    function setPose(i)
        st.qa = poses(i).qa; st.q = poses(i).q;
        set([hPath hRaw], 'XData', nan, 'YData', nan, 'ZData', nan); clearpoints(hTrail);
        applyPose();
    end

    function fromSliders()
        v = arrayfun(@(h) h.Value, sl);
        st.qa = deg2rad(v(1:4)).'; st.q = [deg2rad(v(5:6)).'; v(7)/mm];
        applyPose();
    end

    function planHere()
        if st.busy, return; end
        st.busy = true; cl = onCleanup(@() setBusy(false));
        st.q = P.q_home;
        ttl.String = 'Planning: RRT-Connect in the 4-D arm C-space ...'; drawnow;
        try
            PL = fp_arm_plan(st.qa, qPre, P, struct('dSafe', 25, 'verifyEvery', 4));
        catch err
            ttl.String = ['Planning failed: ' err.message]; return
        end
        S = PL.stats;
        rawF = zeros(3,0);
        for ir = 1:size(PL.raw,2)-1
            for s = linspace(0, 1, 12)
                TF = fp_arm_fk(PL.raw(:,ir) + s*(PL.raw(:,ir+1) - PL.raw(:,ir)), P); rawF(:,end+1) = TF(1:3,4); %#ok<AGROW>
            end
        end
        set(hRaw, 'XData', rawF(1,:)*mm, 'YData', rawF(2,:)*mm, 'ZData', rawF(3,:)*mm);
        set(hPath, 'XData', PL.pF(1,:)*mm, 'YData', PL.pF(2,:)*mm, 'ZData', PL.pF(3,:)*mm);
        clearpoints(hTrail);
        ttl.String = sprintf(['RRT-Connect: %d iterations, %d -> %d waypoints, %.1f s  |  min clearance %.1f mm (%s vs %s)' ...
            '  |  min stop margin %.1f mm'], S.iter, S.nRaw, S.nShort, S.time, S.minClr, S.minBody, S.minObst, S.minSSM);
        for k = 1:size(PL.Q,2)
            if ~isvalid(fig), return; end
            t0 = tic; st.qa = PL.Q(:,k);
            TF = fp_arm_fk(st.qa, P); addpoints(hTrail, TF(1,4)*mm, TF(2,4)*mm, TF(3,4)*mm);
            applyPose(sprintf(' speed %5.1f mm/s  (limit %5.1f)', max(PL.v(k), PL.vR(k))*mm, PL.vlim(k)*mm));
            pause(max(0, 0.04 - toc(t0)));
        end
    end

    function setBusy(b), st.busy = b; end

    function runMain()
        a = questdlg(['The full analysis takes about 15 minutes, closes this window and writes all figures ' ...
            'and results/run_log.txt. Continue?'], 'Run main_final_project', 'Run', 'Cancel', 'Cancel');
        if strcmp(a, 'Run'), evalin('base', 'main_final_project'); end
    end

    function goView(v, dur)
        if nargin < 2, dur = 0.8; end
        d = [sind(v{2})*cosd(v{3}), -cosd(v{2})*cosd(v{3}), sind(v{3})];
        tgt1 = v{1}(:).'; pos1 = tgt1 + v{4}*d; va1 = v{5};
        animateCam(pos1, tgt1, va1, [0 0 1], dur);
    end

    function cameraView()
        K = P.cam; Cc = K.C*mm; tg = K.target*mm; zc = (tg - Cc)/norm(tg - Cc);
        bx = cross(zc, [0;0;1]); bx = bx/norm(bx); by = cross(zc, bx);
        set(H.camHead, 'Visible', 'off');
        animateCam(Cc.' - 2*zc.', tg.', 2*atand((K.res(2)/2)/K.f), -by.', 0.8);
        ttl.String = 'View from the left/right stereo camera mid-point (camera head hidden)';
    end

    function animateCam(pos1, tgt1, va1, up1, dur)
        if ~strcmp(ax.Projection, 'perspective'), ax.Projection = 'perspective'; end
        pos0 = ax.CameraPosition; tgt0 = ax.CameraTarget; va0 = ax.CameraViewAngle; up0 = ax.CameraUpVector;
        if isequal(up1, [0 0 1]) && ~isempty(H.camHead) && isgraphics(H.camHead), set(H.camHead, 'Visible', 'on'); end
        n = max(1, round(dur*30));
        for ic = 1:n
            s = ic/n; s = s*s*(3 - 2*s);
            up = up0 + s*(up1 - up0); if norm(up) < 1e-6, up = up1; end
            set(ax, 'CameraPosition', pos0 + s*(pos1 - pos0), 'CameraTarget', tgt0 + s*(tgt1 - tgt0), ...
                'CameraViewAngle', va0*(va1/va0)^s, 'CameraUpVector', up);
            if dur > 0, drawnow; pause(0.01); end
        end
        drawnow;
    end

    function mouseDown(~, ~)
        if ~isOverAxes(), return; end
        switch fig.SelectionType
            case 'normal', nav.mode = 'orbit';
            case {'alt', 'extend'}, nav.mode = 'pan';
            otherwise, nav.mode = '';
        end
        nav.p0 = get(fig, 'CurrentPoint');
    end

    function mouseUp(~, ~), nav.mode = ''; end

    function mouseMove(~, ~)
        if isempty(nav.mode), return; end
        p = get(fig, 'CurrentPoint'); d = p - nav.p0; nav.p0 = p;
        switch nav.mode
            case 'orbit', camorbit(ax, -0.3*d(1), -0.3*d(2), 'data', [0 0 1]);
            case 'pan',   camdolly(ax, -d(1)/300, -d(2)/300, 0, 'movetarget', 'fixtarget');
        end
        drawnow limitrate;
    end

    function ok = isOverAxes()
        p = get(fig, 'CurrentPoint'); f = getpixelposition(pM, true);
        ok = p(1) >= f(1) && p(1) <= f(1) + f(3) && p(2) >= f(2) && p(2) <= f(2) + f(4);
    end

    function wheel(~, e)
        if ~isOverAxes(), return; end
        camdolly(ax, 0, 0, -0.12*e.VerticalScrollCount, 'fixtarget'); drawnow limitrate;
    end

    function key(~, e)
        f = ax.CameraTarget - ax.CameraPosition; L = norm(f); f = f/L;
        r = cross(f, [0 0 1]); if norm(r) < 1e-6, r = [1 0 0]; end; r = r/norm(r);
        stp = max(20, 0.04*L);
        mv = [0 0 0];
        switch e.Key
            case 'w', mv = stp*f;  case 's', mv = -stp*f;
            case 'a', mv = -stp*r; case 'd', mv = stp*r;
            case 'q', mv = [0 0 -stp]; case 'e', mv = [0 0 stp];
            case 'leftarrow',  campan(ax, 3, 0, 'data', [0 0 1]);
            case 'rightarrow', campan(ax, -3, 0, 'data', [0 0 1]);
            case 'uparrow',    campan(ax, 0, 2, 'data', [0 0 1]);
            case 'downarrow',  campan(ax, 0, -2, 'data', [0 0 1]);
            case 'r', goView(views{2,2});
        end
        if any(mv)
            set(ax, 'CameraPosition', ax.CameraPosition + mv, 'CameraTarget', ax.CameraTarget + mv);
        end
        drawnow limitrate;
    end
end
