function mid = midpoints(mesh)
vert = vertices(mesh);
egs = edges(mesh);
mid = (vert(egs(:,1),:)+vert(egs(:,2),:))/2;
end