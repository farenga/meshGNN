clear, clc, close all

% load the mesh
load('fine2D_rand_voro.mat')

% plot the mesh
plot(mesh)

% get the mesh vertices
V = vertices(mesh);

% get the mesh edges
egs = edges(mesh);

% get the number of elements: [#vertice, #edges, #polygons]
Nvect = mesh.elem_num;

% get the area of the polygons
A = area_faces(mesh);

% get the meshsize of the grid
h = meshsize(mesh);

% create and empty mesh
dim = 2; emtpy_mesh = RMesh(dim);

% Warining: mesh objects are pointers (C style)
% so be careful with assignments
% and consider using copy(mesh) to create a new object

