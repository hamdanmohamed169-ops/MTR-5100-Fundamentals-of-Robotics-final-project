function fp_set_camera(ax, target, az, el, dist, va)
%FP_SET_CAMERA  Place a perspective camera looking at TARGET (mm) from
%   azimuth AZ / elevation EL (deg, same convention as VIEW), at distance
%   DIST (mm) with view angle VA (deg).
d = [sind(az)*cosd(el), -cosd(az)*cosd(el), sind(el)];
set(ax, 'Projection', 'perspective', 'CameraTarget', target(:).', ...
    'CameraPosition', target(:).' + dist*d, 'CameraUpVector', [0 0 1], 'CameraViewAngle', va);
end
