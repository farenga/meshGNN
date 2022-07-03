clear, clc, close all
% GENERATE TRIANGLES DATASET

DatasetSize = 20;
mesh_type = 'tria';

AdjacencyMatrices = {};
AreaVectors = {};
CoordMatrices = {};

for ii=1:DatasetSize
    mesh = generateSortedTriangles(randi([5, 20]));
    mesh_name = strcat(mesh_type,num2str(ii));
    save(['data/', mesh_name],'mesh')

    [graph,W] = connectivity(mesh);
    AdjacencyMatrices{ii} = W>0;
    AreaVectors{ii} = area_faces(mesh)';
    CoordMatrices{ii} = cell2mat(graph.elem(1));

    fprintf(['Mesh ',int2str(ii),' generated\n'])
end

save('data/AdjacencyMatrices_tria.mat','AdjacencyMatrices')
save('data/AreaVectors_tria.mat','AreaVectors')
save('data/CoordMatrices_tria.mat','CoordMatrices')

fprintf(['Dataset generated\n'])