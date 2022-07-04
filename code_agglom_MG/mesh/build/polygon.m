function [mesh] = polygon(N)
pgon = nsidedpoly(N); % center in (0,0) and radius 1
mesh = RMesh(2);
for i = 1:N
    add_new0(mesh,pgon.Vertices(i,:));
end
poly = 1:N;
view = RView(mesh,1:N,0);
add_poly2D(view,poly);
normalize(mesh)
end