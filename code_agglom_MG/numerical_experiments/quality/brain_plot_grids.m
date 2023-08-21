clear, clc, close all
rng('default')

grids = {'brain'};
methods = {'metis','kmeans','GNN'};
names = {'metis','k-means','GNN'};

G = length(grids);
M = length(methods);
matrix = cell(G,M+1);


for g = 1:G
    load(['fine2D_',grids{g}])
    matrix{g,1} = copy(mesh);
end

for g = 1:G
    for m = 1:M
        load([grids{g},'_',methods{m}])
        matrix{g,m+1} = copy(aggl_mesh{end});
    end
end

boom = false;
M = matrix;
matrix = cell(2,2);
matrix(1,1) = M(1,1);
matrix(1,2) = M(1,2);
matrix(2,1) = M(1,3);
matrix(2,2) = M(1,4);
plot_grids(matrix,boom)

%% labels
font = 14;
names = [{'initial'},names];
grids = {''};

for i = 1:4
    nexttile(i)
    title(names{i},'fontweight','bold','fontsize',font)
end


f = gcf;
f.Position = [489 41.8000 748.8000 740.8000];
f.Name = 'uniform_agglomeration_brain';
save_all_figures