clear, clc, close all

% create a voronoi mesh from 10 random points
P = rand(10,2);
mesh = voromesh2D(P);

% compute the area of each polygon
Area = area_faces(mesh);

% get the number of polygons in the mesh
Npoly = mesh.elem_num(end);

Perim = zeros(1,Npoly);
for i = 1:Npoly
    poly = RView(mesh,i,2); % for each polygon
    Perim(i) = sum(length_edges(poly)); % sum the length of its edges
end
