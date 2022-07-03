function [polygon,vert_ids] = edges2polygon(vertices,edges)
n_edges = size(edges,1);
vert_ids = zeros(n_edges,1);

vert_ids(1) = edges(1,1);
v = edges(1,2);
edges(1,:) = 0;
for i = 2:n_edges
    vert_ids(i) = v;
    [row,col] = find(v == edges);
    new_col = mod(col,2)+1; % cyclic index [1 2]
    new_v = edges(row,new_col);
    edges(row,:) = 0;
    v = new_v;
end
polygon = vertices(vert_ids,:);
end