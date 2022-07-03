clear, clc

grids = {'tria','rand_tria','rand_voro','quads'};
methods = {'metis','kmeans','GNN_base','GNN_Res'};

G = length(grids);
M = length(methods);
iter = zeros(G,M);

for p = 2:3
    for smooth = 5:10
        fprintf(['p = ',num2str(p),', smooth = ',num2str(smooth),'\n'])
        for g = 1:G
            for m = 1:M
                name = [grids{g},'_',methods{m}];
                disp(name)
                load(name,'aggl_mesh');
                rng('default')
                iter(g,m) = multigrid_iterations(aggl_mesh,smooth,p);
                fprintf('\n')
            end
        end
        save([path2('MG'),'MG_iterations_base_p',num2str(p),'_s',num2str(smooth),'_m',num2str(m)],'methods','grids','iter')
        
        T = table(iter(:,1),iter(:,2),iter(:,3),iter(:,4),'VariableNames',methods,'RowNames',grids);
        
        disp(T)
    end
end