clc, clear, close all
load("brain_no_holes.mat")

figure('Name','brain_original')
plot(mesh)
aggl_mesh_brain{1} = mesh;

% Use the kmeans algorithm to perform agglomeration:

% agglomeration strategy (can be a neural network)
% aggl_fun = @aggl_kmeans_fun;
% method = 'kmeans';
aggl_fun = @aggl_GNN_fun;
method = 'GNN';
% aggl_fun = @aggl_metis_fun;
% method = 'metis';

% vector of different target sizes for the agglomerated grids

h_list = max(pdist(vertices(mesh)))*[0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2]; % in descending order!
N_grids = length(h_list);

% compute the indexes of the polygons to agglomerate
new_elem_list = agglom_strategy_vect(mesh,h_list,aggl_fun);

%% compute the agglomerated grids and plot them
for i = 1:N_grids
    aggl_mesh = agglomerate(mesh,new_elem_list{i});

    figure('Name',['brain_agglom_',method,num2str(i)]), plot(aggl_mesh)
    aggl_mesh_brain{i+1} = aggl_mesh;
    % title(['h = ',num2str(h_list(i))])
end

save([path3('aggl_mesh_brain'),'aggl_mesh_brain_',method,'base'],'aggl_mesh_brain')

save_all_figures