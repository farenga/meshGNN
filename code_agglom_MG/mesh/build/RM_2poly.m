function [ids,loc_vert] = RM_2poly(polyview)
egs = edges(polyview);
n_edges = size(egs,1);
ids = zeros(1,n_edges);
ids(1) = egs(1,1);
v = egs(1,2);
egs(1,:) = 0;
for i = 2:n_edges
    ids(i) = v;
    [row,col] = find(v == egs);
    new_col = mod(col,2)+1; % cyclic index [1 2]
    new_v = egs(row,new_col);
    egs(row,:) = 0;
    v = new_v;
end
loc_vert = local_basis(polyview);
if not(ispolycw(loc_vert(ids,1),loc_vert(ids,2)))
	ids = ids(end:-1:1);
end
loc_vert = loc_vert(ids,:);
end