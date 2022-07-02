function L = length_edges(mesh)
vert = vertices(mesh);
egs = edges(mesh);
L = vecnorm(vert(egs(:,1),:)-vert(egs(:,2),:),2,2);
end