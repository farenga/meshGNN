function IDX = fix_partition(mesh,IDX,weights)


mesh1 = RView(mesh,find(IDX==1),mesh.dim);
mesh2 = RView(mesh,find(IDX==2),mesh.dim);
graph1 = connectivity(mesh1,mesh1.dim);
graph2 = connectivity(mesh2,mesh2.dim);
[clust1,n_clust1] = clusters(graph1);
[clust2,n_clust2] = clusters(graph2);

if n_clust1 == 1 && n_clust2 == 1
    return
end

W1 = zeros(1,n_clust1);
for i = 1:n_clust1
    W1(i) = sum(weights(mesh1.loc2glob{mesh.dim+1}(clust1==i)));
end

W2 = zeros(1,n_clust2);
for i = 1:n_clust2
    W2(i) = sum(weights(mesh2.loc2glob{mesh.dim+1}(clust2==i)));
end
 
[m1,I1]=min(W1);
[m2,I2]=min(W2);

if n_clust1 == 1 || (m1 > m2 && n_clust2 > 1)
    IDX(mesh2.loc2glob{mesh.dim+1}(clust2==I2)) = 1;
    IDX = fix_partition(mesh,IDX,weights);
    return
end

if n_clust2 == 1 || (m2 > m1 && n_clust1 > 1)
    IDX(mesh1.loc2glob{mesh.dim+1}(clust1==I1)) = 2;
    IDX = fix_partition(mesh,IDX,weights);
    return
end


% if M1 > M2
%     IDX(:) = 2;
%     IDX(mesh1.loc2glob{mesh.dim+1}(clust1==I1)) = 1;
% else
%     IDX(:) = 1;
%     IDX(mesh2.loc2glob{mesh.dim+1}(clust2==I2)) = 2;
% end


% egs = edges(graph);
% 
% 
% % problema ci sono 2 cluster ma isolati <----
% IDX_old = IDX;
% 
% if nnz(IDX==1) > 1 && nnz(IDX==2) > 1
%     ok = 0;
%     while ok == 0
%         ok = 1;
%         for i = 1:graph.elem_num(1)
%             neigh =[egs(egs(:,1)==i,2);egs(egs(:,2)==i,1)];
% 
%             if all(IDX(neigh)~=IDX(i))
%                 ok = 0;
%                 if IDX(i) == 1
%                     IDX(i) = 2;
%                 else
%                     IDX(i) = 1;
%                 end
%             end
%         end
%     end
% end
% 
% if all(IDX==IDX(1))
%     I_alone = find(IDX_old~=IDX(1));
%     [~,IWmax] = max(weights(I_alone));
%     IDX(I_alone(IWmax)) = IDX_old(I_alone(IWmax));
% end