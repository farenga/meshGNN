clear, clc, close all

% create a voronoi mesh from 10 random points
P = rand(10,2);
mesh = voromesh2D(P);

% plot the mesh and show the number of each polygon
% polygon   -> dim_numbers = 2
% egdes     -> dim_numbers = 1
% vertices  -> dim_numbers = 0
dim_numbers = 2;
plot(mesh,dim_numbers)

% compute the connecitivity graph:
% each polygon is a node
% each edges connects two neighbouring polygons
[graph,W] = connectivity(mesh);

% W is the adjacency matrix:
% W(i,j) = 0 if there is NO connection between polygons i and j
% W(i,j) > 0 if there is connection between polygons i and j
% W(i,j) is the length of the common edges between polygons i and j

% graph has the same structure of a mesh object
% plot the graph on top of the mesh
figure, plot(mesh)
hold on, plot(graph)

