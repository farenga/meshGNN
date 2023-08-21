clear, clc

grids = {'tria','rand_tria','rand_voro','quads'};
methods = {'metis','kmeans','GNN'};

G = length(grids);
M = length(methods);
iter = zeros(G,M,2);

for L = 2 % list=[2,3,4], default=3
    for p = 1 % list=[1,2,3], default=1
        for m = 5 % list=[1,3,5], default=3
            filename = ['MG_iter_p',num2str(p),'_m',num2str(m),'_L',num2str(L),'.mat'];
            if not(exist(filename,'file'))
            fprintf(['p = ',num2str(p),', m = ',num2str(m),', L = ',num2str(L),'\n'])
            for i = 1:G
                for j = 1:M
                    name = [grids{i},'_',methods{j}];
                    disp(name)
                    load(name,'aggl_mesh');
                    rng('default')
                    [iter(i,j,1),iter(i,j,2)] = multigrid_iterations(aggl_mesh(1:L),m,p);
                    fprintf('\n')
                end
            end
            save([path2('MG'),filename],'methods','grids','iter')
            end
            
            T = table(iter(:,1),iter(:,2),iter(:,3),'VariableNames',methods,'RowNames',grids);
            
            disp(T)
        end
    end
end