clear, clc, close all

methods = {'metis','kmeans','GNN_base','GNN_Res'};

M = length(methods);

for j = 1:M
    fprintf(['brain',' ',methods{j},'\n'])
    load([path3('aggl_mesh_brain'),'aggl_mesh_brain_',methods{j},'.mat'])
    [UF,CR] = quality(aggl_mesh_brain{end});
    save([path3('aggl_mesh_brain'),'quality_brain_',methods{j}],'UF','CR');
end 
