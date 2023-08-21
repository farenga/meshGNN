clear, clc, close all
rng('default')


h_list = sqrt(2)./[2 4 8];
n_aggl = length(h_list);

grids = {'rand_voro'};
% grids = {'tria','rand_tria','rand_voro','quads'};
% methods = {'kmeans','GNN_base','GNN_res'}; % Attention: run also GNN_Res, by changing the aggl_mesh_GNN function

methods = {'kmeans','GNN'};

% net = 'res';
% methods = {['GNN_',net]};


for g = 1:length(grids)
    for m = 1:length(methods)
        display([grids{g},' ',methods{m}])
        load([path2('grids'),'fine2D_',grids{g}])
        
%         h0 = meshsize(mesh);
%         %%%
%         h_list = [h0*4,h0*2];
%         n_aggl = length(h_list);
%         %%%

        aggl_mesh = cell(1,n_aggl+1);
        aggl_mesh{1} = mesh;
      
        tic
        switch methods{m}
            case 'metis'
                error('run metis_vanilla.m')
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

        for i = 1:n_aggl+1
            H(i) = meshsize(aggl_mesh{i});
        end

        save([path2('grids'),grids{g},'_',methods{m}],...
            'aggl_mesh','aggl_time')
    end
end

