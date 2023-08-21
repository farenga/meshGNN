clear, clc, close all
rng('default')

load("fine2D_brain.mat");
D = max(pdist(vertices(mesh)));
h_list = 0.2*D;
n_aggl = length(h_list);

grids = {'brain'};
methods = {'kmeans','GNN'};


for g = 1:length(grids)
    for m = 1:length(methods)
        display([grids{g},' ',methods{m}])
        load([path2('grids'),'fine2D_',grids{g}])
        
        aggl_mesh = cell(1,n_aggl+1);
        aggl_mesh{1} = mesh;
      
        tic
        switch methods{m}
            case 'metis'
                error('run vanilla_metis.m')
                %aggl_fun = @aggl_metis_fun;
            case 'kmeans'
                aggl_fun = @aggl_kmeans_fun;
            case 'GNN'
                aggl_fun = @aggl_GNN_fun;
        end

        new_elem_list= agglom_strategy_vect(mesh,h_list,aggl_fun);
        aggl_time = toc;
        
        for i = 1:n_aggl
            aggl_mesh{n_aggl+1-i+1} = agglomerate(mesh,new_elem_list{i});
        end

        save([path2('grids'),grids{g},'_',methods{m}],...
            'aggl_mesh','aggl_time')
    end
end

