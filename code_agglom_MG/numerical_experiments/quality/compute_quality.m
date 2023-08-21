clear, clc, close all

%grids = {'tria','rand_tria','rand_voro','quads'}; % aggl_mesh{end-1} in line 13
grids = {'brain'}; % aggl_mesh{end} in line 13
methods = {'metis','kmeans','GNN'};

M = length(methods);
G = length(grids);
for i = 1:G
    for j = 1:M
        fprintf([grids{i},' ',methods{j},'\n'])
        load([path2('grids'),grids{i},'_',methods{j},'.mat'])
        [UF,CR] = quality(aggl_mesh{end});
%         UF = []; CR = [];
%         for k = 1:length(aggl_mesh)
%         [UF_k,CR_k] = quality(aggl_mesh{k});
%         UF = [UF,UF_k];
%         CR = [CR,CR_k];
%         end

        save([path2('grids'),'quality_',grids{i},'_',methods{j}],'UF','CR');
    end 
end
