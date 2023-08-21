clear, clc, close all
rng('default')

grids = {'tria','rand_tria','rand_voro','quads'};
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
        matrix{g,m+1} = copy(aggl_mesh{end-1});
    end
end

boom = false;
plot_grids(matrix,boom)

%% labels
font = 14;
names = [{'initial'},names];
grids = {'triangles','random','voronoi','squares'};

for m = 1:M+1
    nexttile(m)
    title(names{m},'fontweight','bold','fontsize',font)
end
for g = 1:G
    nexttile((g-1)*(M+1)+1);
    ax = gca;
    ax.YLabel.Visible = 'on';
    ylabel(grids{g},'fontweight','bold','fontsize',font)
end

f = gcf;
f.Position = [489 41.8000 748.8000 740.8000];
f.Name = 'uniform_agglomeration';
save_all_figures