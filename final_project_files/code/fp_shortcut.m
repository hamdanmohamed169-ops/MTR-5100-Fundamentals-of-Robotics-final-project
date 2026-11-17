function path = fp_shortcut(isFree, path, nTry, res)
%FP_SHORTCUT  Path smoothing by random shortcuts (Module 8):
%   pick two waypoints, try the straight segment, keep it only if every
%   intermediate configuration is collision-free.  path is 2xK.
for k = 1:nTry
    K = size(path,2);
    if K <= 2, break; end
    ij = sort(randperm(K, 2));
    i = ij(1); j = ij(2);
    if j - i < 2, continue; end
    a = path(:,i); b = path(:,j);
    n = max(2, ceil(norm(b-a)/res) + 1);
    ok = true;
    for s = linspace(0,1,n)
        if ~isFree(a + s*(b-a)), ok = false; break; end
    end
    if ok, path = [path(:,1:i) path(:,j:end)]; end
end
end
