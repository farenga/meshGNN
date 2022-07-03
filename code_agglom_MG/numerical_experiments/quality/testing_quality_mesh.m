clear, clc, close all

grids = {'tria','randtria','voro','quad'};
methods = {'metis','kmeans','GNN_base','GNN_Res'};

N = 50;

M = length(methods);
G = length(grids);

quality_UF = {zeros(N,length(grids)),zeros(N,length(grids)),zeros(N,length(grids))};
quality_CR = {zeros(N,length(grids)),zeros(N,length(grids)),zeros(N,length(grids))};

for ii = 1:N
    fprintf(['N = ',num2str(ii),'\n'])
    for i = 1:G
        for j = 1:M
            fprintf([grids{i},' ',methods{j},'\n'])
            load([path2('testing'),grids{i},'_',methods{j},num2str(ii),'.mat'])
            [UF,CR] = quality(aggl_mesh{end});
            quality_UF{j}(ii,i) = mean(UF);
            quality_CR{j}(ii,i) = mean(CR);
        end 
    end
end

% Compute average
avg_vector_UF = zeros(G,M);
avg_vector_CR = zeros(G,M);
for i = 1:G
    for j = 1:M
        avg_vector_UF(i,j) = sum(quality_UF{j}(:,i))/N;
        avg_vector_CR(i,j) = sum(quality_CR{j}(:,i))/N;
    end 
end

save([path2('testing'),'quality_averages'],'avg_vector_UF','avg_vector_CR');