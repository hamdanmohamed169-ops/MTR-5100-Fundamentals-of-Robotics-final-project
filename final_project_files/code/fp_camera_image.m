function [img, pp, hw] = fp_camera_image(P, REG, qa)
%FP_CAMERA_IMAGE  Synthetic frame of the LEFT tracking camera.
%   Renders the 3-D prototype scene (FP_SCENE) from the left camera's exact
%   pose REG.cams(1) with the pinhole intrinsics of P.cam (focal length f px,
%   principal point at the image centre). The arm is placed at joints qa.
%   Returns an RGB image img covering pixels pp +/- hw (u and v), so the
%   detected limbus points REG.uv{1} and the fits overlay in the same pixel
%   coordinates:  image(ax, 'XData', pp(1)+[-hw hw], 'YData', pp(2)+[-hw hw], 'CData', img)
hw = 185;                                          % half window (px) around the principal point
pp = P.cam.res/2;
fig2 = figure('Visible', 'off', 'Color', 'k', 'Position', [40 40 720 720]);
ax2 = axes(fig2, 'Position', [0 0 1 1], 'Color', 'k');
H2 = fp_scene(ax2, P);
fp_scene_update(H2, P.q_home, P, 0, 0, qa);
set([H2.pEmark H2.camHead], 'Visible', 'off');      % no overlay marker; we look out from inside the camera head
axis(ax2, 'off'); daspect(ax2, [1 1 1]);
xlim(ax2, [-45 45]); ylim(ax2, [-45 45]); zlim(ax2, [-25 60]);   % tight box -> fine depth-buffer precision
cam = REG.cams(1); Cc = cam.Cc.'*1e3;
light(ax2, 'Position', [0.4 -1 2], 'Style', 'infinite', 'Color', [1 0.97 0.92]);
light(ax2, 'Position', [-1 0.8 0.5], 'Style', 'infinite', 'Color', [0.45 0.5 0.6]);
lighting(ax2, 'gouraud'); ax2.AmbientLightColor = [0.85 0.88 0.95];
% camera y-axis points down the image, so "up" on screen is -y_cam
set(ax2, 'Projection', 'perspective', 'CameraPosition', Cc, 'CameraTarget', P.cam.target.'*1e3, ...
    'CameraUpVector', -cam.R(:,2).', 'CameraViewAngle', 2*atand(hw/P.cam.f));
drawnow;
Fr = getframe(ax2); img = double(Fr.cdata)/255;
close(fig2);
[n1, n2, ~] = size(img);                           % mild lens vignetting + sensor noise
[gx, gy] = meshgrid(linspace(-1, 1, n2), linspace(-1, 1, n1));
s = RandStream('twister', 'Seed', 3);
img = min(max(img.*(1 - 0.18*(gx.^2 + gy.^2)) + 0.012*randn(s, n1, n2), 0), 1);
end
