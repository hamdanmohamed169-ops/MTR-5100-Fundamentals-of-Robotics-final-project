function [path, info] = fp_rrt_connect(isFree, qs, qg, lo, hi, opt)
%FP_RRT_CONNECT  Bidirectional RRT (RRT-Connect, Kuffner & LaValle 2000) in an
%   n-D configuration space (Module 8).
%   isFree(q)  : collision test of one configuration (n x 1)
%   qs, qg     : start / goal (n x 1);  lo, hi : joint limits (n x 1)
%   opt.step   : extension step (rad);   opt.edgeRes : edge check resolution (rad)
%   opt.maxIter: iterations;  opt.w : per-joint metric weights (n x 1)
%   One tree grows from the start and one from the goal; each iteration
%   extends one tree towards a random sample and then greedily CONNECTs the
%   other tree to the new node; the trees swap roles every iteration.
%   Every accepted edge is collision-checked every opt.edgeRes, or by
%   opt.edgeFree(a, b) -> [ok, nChecks] when it is given.
n = numel(qs); w = opt.w(:);
dist = @(A, b) sqrt(sum((w.*(A - b)).^2, 1));
Ta = struct('V', qs(:), 'par', 0); Tb = struct('V', qg(:), 'par', 0);
info = struct('iter', 0, 'success', false, 'nodes', 0, 'checks', 0);
path = [];
if ~isFree(qs) || ~isFree(qg), return; end
fromStart = true;
for it = 1:opt.maxIter
    info.iter = it;
    qr = lo + (hi - lo).*rand(n, 1);
    [Ta, ka, st] = extend(Ta, qr);
    if st ~= 0                                          % new node added -> try to connect the other tree
        [Tb, kb, st2] = connect(Tb, Ta.V(:,ka));
        if st2 == 2
            pa = branch(Ta, ka); pb = branch(Tb, kb);
            if fromStart, path = [pa fliplr(pb(:,1:end-1))]; else, path = [pb fliplr(pa(:,1:end-1))]; end
            path = path(:, [true any(diff(path, 1, 2) ~= 0, 1)]);
            info.success = true; break
        end
    end
    [Ta, Tb] = deal(Tb, Ta); fromStart = ~fromStart;     % swap the trees
end
info.nodes = size(Ta.V,2) + size(Tb.V,2);
if info.success
    if any(abs(path(:,1) - qs(:)) > 1e-12), path = fliplr(path); end
end

    function [T, k, status] = extend(T, q)
        % status 0 trapped, 1 advanced, 2 reached q
        [~, j] = min(dist(T.V, q)); qn = T.V(:,j);
        dq = q - qn; L = sqrt(sum((w.*dq).^2));
        if L < 1e-12, k = j; status = 2; return; end
        if L <= opt.step, qnew = q; status = 2; else, qnew = qn + dq*(opt.step/L); status = 1; end
        if ~edgeFree(qn, qnew), k = j; status = 0; return; end
        T.V(:,end+1) = qnew; T.par(end+1) = j; k = size(T.V,2);
    end

    function [T, k, status] = connect(T, q)
        status = 1;
        while status == 1
            [T, k, status] = extend(T, q);
        end
    end

    function ok = edgeFree(a, b)
        if isfield(opt, 'edgeFree') && ~isempty(opt.edgeFree)      % caller-supplied (e.g. adaptive) edge check
            [ok, nc] = opt.edgeFree(a, b); info.checks = info.checks + nc; return
        end
        m = max(2, ceil(max(abs(b - a))/opt.edgeRes) + 1);
        ok = true;
        for s = linspace(0, 1, m)
            info.checks = info.checks + 1;
            if ~isFree(a + s*(b - a)), ok = false; return; end
        end
    end
end

function p = branch(T, k)
% nodes from the root of tree T to node k
p = T.V(:,k);
while T.par(k) ~= 0
    k = T.par(k); p = [T.V(:,k) p]; %#ok<AGROW>
end
end
