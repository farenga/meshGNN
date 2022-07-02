clear, clc, close all

% create a voronoi mesh from 10 random points
P = rand(20,2);
mesh = voromesh2D(P);
%mesh = load('saved_mesh/mesh_100.mat').mesh;
figure, plot(mesh), title('original mesh')

% Use the kmeans algorithm to perform agglomeration:

% agglomeration strategy (can be a neural network)
%aggl_fun = @aggl_kmeans_fun;
aggl_fun = @aggl_GNN_fun;

% vector of different target sizes for the agglomerated grids
h_list = [0.9,0.7,0.5]; % in descending order!
N_grids = length(h_list);

% compute the indexes of the polygons to agglomerate
new_elem_list = agglom_strategy_vect(mesh,h_list,aggl_fun);

% compute the agglomerated grids and plot them
for i = 1:N_grids
    aggl_mesh = agglomerate(mesh,new_elem_list{i});
    figure, plot(aggl_mesh)
    title(['h = ',num2str(h_list(i))])
end