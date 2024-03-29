clear, clc, close all
rng('default')

%%%
h_list = sqrt(2)./[4 8];
n_aggl = length(h_list);
%%%

N = 1;
grids = {'tria','randtria','voro','quad'};
methods = {'GNN'}; % Attention: run also GNN_Res, by changing the aggl_mesh_GNN function

for ii = 1:50
    display(['N = ',num2str(ii)])
    for g = 1:length(grids)
        for m = 1:length(methods)
            display([grids{g},' ',methods{m}])
            load([path2('testing'),grids{g},num2str(ii)])
            aggl_mesh = cell(1,n_aggl+1);
            aggl_mesh{1} = mesh;
            tic
            switch methods{m}
                case 'metis'
                    aggl_fun = @aggl_metis_fun;
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
    
            save([path2('testing'),grids{g},'_',methods{m},num2str(ii)],...
                'aggl_mesh','aggl_time')
        end
    end
end