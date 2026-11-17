function [path, info] = fp_astar(occ, s, g)
%FP_ASTAR  A* search on an 8-connected C-space grid (Module 8).
%   occ  : logical grid, true = configuration in collision (rows = q2, cols = q1)
%   s, g : [row col] start / goal cells
%   g-cost = travelled Euclidean grid distance; h = Euclidean distance to goal
%   (admissible and consistent -> the returned path is optimal on the grid).
[nr, nc] = size(occ);
idx = @(r,c) (c-1)*nr + r;
G = inf(nr,nc); F = inf(nr,nc); parent = zeros(nr,nc);
closed = false(nr,nc);
G(s(1),s(2)) = 0; F(s(1),s(2)) = hypot(s(1)-g(1), s(2)-g(2));
openList = false(nr,nc); openList(s(1),s(2)) = true;
moves = [-1 -1; -1 0; -1 1; 0 -1; 0 1; 1 -1; 1 0; 1 1];
info.expanded = 0;
while any(openList(:))
    Fo = F; Fo(~openList) = inf;
    [~, k] = min(Fo(:));
    [r, c] = ind2sub([nr nc], k);
    if r == g(1) && c == g(2), break; end
    openList(r,c) = false; closed(r,c) = true;
    info.expanded = info.expanded + 1;
    for m = 1:8
        rr = r + moves(m,1); cc = c + moves(m,2);
        if rr < 1 || rr > nr || cc < 1 || cc > nc, continue; end
        if occ(rr,cc) || closed(rr,cc), continue; end
        % forbid diagonal corner cutting through obstacles
        if moves(m,1)~=0 && moves(m,2)~=0 && (occ(r,cc) || occ(rr,c)), continue; end
        gn = G(r,c) + norm(moves(m,:));
        if gn < G(rr,cc)
            G(rr,cc) = gn; F(rr,cc) = gn + hypot(rr-g(1), cc-g(2));
            parent(rr,cc) = idx(r,c); openList(rr,cc) = true;
        end
    end
end
if isinf(G(g(1),g(2)))
    path = []; info.cost = inf; return
end
path = g; k = parent(g(1),g(2));
while k ~= 0
    [r,c] = ind2sub([nr nc], k);
    path = [r c; path]; %#ok<AGROW>
    k = parent(r,c);
end
info.cost = G(g(1),g(2));
end
