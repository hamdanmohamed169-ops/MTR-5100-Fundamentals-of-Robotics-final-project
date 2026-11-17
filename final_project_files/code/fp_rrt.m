function [path, info] = fp_rrt(isFree, qs, qg, lo, hi, opt)
%FP_RRT  Basic RRT in a 2-D C-space (Module 8).
%   isFree(q) : collision test of a single configuration (column 2x1)
%   qs, qg    : start / goal (2x1); lo, hi: C-space bounds
%   opt.step, opt.goalBias, opt.maxIter, opt.edgeRes, opt.goalTol
%   Each proposed extension is accepted only if the whole edge is free
%   (checked every opt.edgeRes).
nodes = qs(:); par = 0; path = [];
info.iter = 0; info.success = false;
for it = 1:opt.maxIter
    info.iter = it;
    if rand < opt.goalBias, qr = qg; else, qr = lo + (hi-lo).*rand(2,1); end
    [~, j] = min(sum((nodes - qr).^2, 1));
    qn = nodes(:,j); d = qr - qn; L = norm(d);
    if L < 1e-12, continue; end
    qnew = qn + d/L*min(opt.step, L);
    if ~edgeFree(isFree, qn, qnew, opt.edgeRes), continue; end
    nodes(:,end+1) = qnew; par(end+1) = j; %#ok<AGROW>
    if norm(qnew - qg) <= opt.goalTol && edgeFree(isFree, qnew, qg, opt.edgeRes)
        nodes(:,end+1) = qg; par(end+1) = size(nodes,2)-1; %#ok<AGROW>
        info.success = true; break
    end
end
info.nodes = nodes; info.parent = par;
if info.success
    k = size(nodes,2);
    while k ~= 0
        path = [nodes(:,k) path]; %#ok<AGROW>
        k = par(k);
    end
end
end

function ok = edgeFree(isFree, a, b, res)
n = max(2, ceil(norm(b-a)/res) + 1);
ok = true;
for s = linspace(0,1,n)
    if ~isFree(a + s*(b-a)), ok = false; return; end
end
end
