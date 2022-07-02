function IDX = aggl_kmeans_fun(mesh,weights)

C = zeros(mesh.elem_num(end),mesh.dim);
for i = 1:mesh.elem_num(end)
    elem = RView(mesh,i,mesh.dim);
    if mesh.dim == 2
        polygon = edges2polygon(vertices(elem),edges(elem));
        P = polyshape(polygon);
        [C(i,1),C(i,2)] = centroid(P);
    end
    if mesh.dim == 3
        V = vertices(elem);
        C(i,:) = conv_centr(V);
    end
end

scores = min(10,round(weights/min(weights)));
C_weight = zeros(sum(scores),mesh.dim);
j = 1; % first empty position
for i = 1:mesh.elem_num(end)
   C_weight(j:(j+scores(i)-1),:)= repmat(C(i,:),[scores(i),1]);
   j = j+scores(i);
end

n_clust = 2;
expanded_IDX = kmeans(C_weight,n_clust,...
    'Replicates',10,...
    'OnlinePhase','on',...
    'Distance','sqeuclidean',...
    'MaxIter',100);

IDX = zeros(1,mesh.elem_num(end));
I = cumsum(scores);
for i = 1:mesh.elem_num(end)
    IDX(i)= expanded_IDX(I(i));
end
