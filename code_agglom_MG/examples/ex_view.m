clear, clc, close all

% load the mesh
load('fine2D_rand_voro.mat')

% plot the mesh
figure, plot(mesh)

% A view allows to work with only some elements of the mesh,
% while hiding the other ones.

% create a view with only polygons from 10 to 20
ids = 10:40; dim = 2;
view = RView(mesh,ids,dim);

% plot the view
figure, plot(view)

% get the nunmber of polygons in the view
N = view.elem_num(3);

% create a new mesh from the view
new_mesh = copy(view);

% get the vertices of polygon 4 of the mesh
id = 4; dim = 2;
V4 = vertices(RView(mesh,id,dim));

% Warining: mesh and view objects are pointers (C style)
% so be careful with assignments
% and consider using copy to create a new object
