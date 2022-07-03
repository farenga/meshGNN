clear, clc, close all
% GENERATE QUAD DATASET

DatasetSize = 20;
mesh_type = 'quad';

AdjacencyMatrices = {};
AreaVectors = {};
CoordMatrices = {};

for ii=1:DatasetSize
    n = randi([3,6],1);
    mesh = generateQuads(n);
    mesh_name = strcat(mesh_type,num2str(ii));
    save(['data/', mesh_name],'mesh')

    [graph,W] = connectivity(mesh);
    AdjacencyMatrices{ii} = W>0;
    AreaVectors{ii} = area_faces(mesh)';
    CoordMatrices{ii} = cell2mat(graph.elem(1));

    fprintf(['Mesh ',int2str(ii),' generated\n'])
end

save('data/AdjacencyMatrices_quad.mat','AdjacencyMatrices')
save('data/AreaVectors_quad.mat','AreaVectors')
save('data/CoordMatrices_quad.mat','CoordMatrices')

fprintf(['Dataset generated\n'])