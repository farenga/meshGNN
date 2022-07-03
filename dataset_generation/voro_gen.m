clear, clc, close all
% GENERATE RANDOM VORONOI DATASET

DatasetSize = 20;
mesh_type = 'rand_voro';

AdjacencyMatrices = {};
AreaVectors = {};
CoordMatrices = {};

for ii=1:DatasetSize
    P = rand(randi([50,1000],1),2);
    mesh = voromesh2D(P);
    mesh_name = strcat(mesh_type,num2str(ii));
    save(['data/', mesh_name],'mesh')

    [graph,W] = connectivity(mesh);
    AdjacencyMatrices{ii} = W>0;
    AreaVectors{ii} = area_faces(mesh)';
    CoordMatrices{ii} = cell2mat(graph.elem(1));

    fprintf(['Mesh ',int2str(ii),' generated\n'])
end

save('data/AdjacencyMatrices_rand_voro.mat','AdjacencyMatrices')
save('data/AreaVectors_rand_voro.mat','AreaVectors')
save('data/CoordMatrices_rand_voro.mat','CoordMatrices')

fprintf(['Dataset generated\n'])