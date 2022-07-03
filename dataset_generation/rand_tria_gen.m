clear, clc, close all
% GENERATE RANDOM TRIANGLES DATASET

DatasetSize = 20;
mesh_type = 'rand_tria';

AdjacencyMatrices = {};
AreaVectors = {};
CoordMatrices = {};

for ii=1:DatasetSize
    mesh = generateRandomTriangles(randi([10,400]));
    mesh_name = strcat(mesh_type,num2str(ii));
    save(['data/', mesh_name],'mesh')

    [graph,W] = connectivity(mesh);
    AdjacencyMatrices{ii} = W>0;
    AreaVectors{ii} = area_faces(mesh)';
    CoordMatrices{ii} = cell2mat(graph.elem(1));

    fprintf(['Mesh ',int2str(ii),' generated\n'])
end

save('data/AdjacencyMatrices_rand_tria.mat','AdjacencyMatrices')
save('data/AreaVectors_rand_tria.mat','AreaVectors')
save('data/CoordMatrices_rand_tria.mat','CoordMatrices')

fprintf(['Dataset generated\n'])