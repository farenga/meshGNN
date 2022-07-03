clear, clc, close all

grids = {'tria','rand_tria','rand_voro','quads'};
methods = {'metis','kmeans','GNN_base','GNN_Res','GNN_base_voro','GNN_Res_voro'};

M = length(methods);
G = length(grids);
for i = 1:G
    for j = 1:M
        fprintf([grids{i},' ',methods{j},'\n'])
        load([path2('grids'),grids{i},'_',methods{j},'.mat'])
        [UF,CR] = quality(aggl_mesh{end});
        save([path2('grids'),'quality_',grids{i},'_',methods{j}],'UF','CR');
    end 
end
