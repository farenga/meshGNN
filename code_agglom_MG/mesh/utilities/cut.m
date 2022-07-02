function cut(mesh,n,x0)
n = n/norm(n);
% refine edges

% n*(x-x0) = 0
% clust1 -n-> clust2
dtol = mesh.tol*mesh.dim;
vert = vertices(mesh);
clust = zeros(1,mesh.elem_num(1));
clust(n*(vert-x0)' < -dtol) = 1;
clust(n*(vert-x0)' > +dtol) = 2;

% if mesh is already on one side of the plane just return
if not(all(ismember([1,2],clust)))
    return
end

clust(clust == 0) = 2;

egs = edges(mesh);
ids = find(clust(egs(:,1))~=clust(egs(:,2)));
v1 = get0(mesh,egs(ids,1));
v2 = get0(mesh,egs(ids,2));

% plane: n*(x-x0) = 0
% line: x = (v2-v1)*t + v1
% dot((v2-v1)*t+v1-x0,n) = 0 --> t = dot(x0-v1,n)/dot(v2-v1,n)
t = ((x0-v1)*n') ./ ((v2-v1)*n');
c = (v2-v1).*t + v1;

on_v1 = vecnorm(c-v1,inf,2)<mesh.tol;
on_v2 = vecnorm(c-v2,inf,2)<mesh.tol;
new_points =  find(not(on_v1) & not(on_v2))';

if min(pdist(c(new_points,:))) == 0
    error('two new points are equal: ill-shaped element')
end

for i = new_points
    ci = add_new0(mesh,c(i,:));
%     ci = add0(RView(mesh),c(i,:));
    ei = egs(ids(i),:);
    add_fun = @(j) add_new1(mesh,[ei(clust(ei)==j),ci]);
    replace1(mesh,ids(i),add_fun,2);
end

new_clusters = 2*ones(1,mesh.elem_num(2));
new_clusters(any(clust(egs)==1,2)) = 1;
new_clusters(ids(on_v1)) = clust(egs(ids(on_v1),2));
new_clusters(ids(on_v2)) = clust(egs(ids(on_v2),1));

% to perform a check:
% plot(connectivity(mesh,1),new_clusters)

% refine upper levels

for dim = 2:manifold(mesh)
    clust = new_clusters;    
    % allocate memory by excess
    new_clusters = zeros(1,2*mesh.elem_num(dim+1));     
    for i = 1:mesh.elem_num(dim+1)
        view1 = RView(mesh,i,dim);
        clust1 = clust(view1.loc2glob{dim});
        ids_clust = unique(clust1);
        if length(ids_clust) > 1
            ids_clust = cut_fun(view1,clust1);
            n_clust = length(ids_clust);
            new_ids = (mesh.elem_num(dim+1)-n_clust+2):mesh.elem_num(dim+1);
            new_clusters(new_ids) = ids_clust(2:end);
        end
        new_clusters(i) =  ids_clust(1);
    end
end

end

function ids_clust = cut_fun(view1,clust1)
% view1 must contain only ONE element

    dim = manifold(view1);
    ids_clust = unique(clust1);
    dtol = view1.tol*view1.dim;
    
    part = cell(1,2);    
    degenerate = zeros(1,2);
    elem = @(j) find(clust1 == ids_clust(j));
    for j = 1:2
        part{j} = RView(view1,elem(j),dim-1);
        vert = vertices(part{j});
        degenerate(j) = rank(vert(2:end,:)-vert(1,:),dtol) < dim;
    end
    if any(degenerate)
        ids_clust = ids_clust(not(degenerate)); 
        return
    end
    
    % add the separating plane
    bp1 = boundary(part{1});
    if dim > 2
        add_new2(bp1,1:bp1.elem_num(dim-1),dim-1);
    else
        add_new1(bp1,1:bp1.elem_num(dim-1));
    end
    
    add_fun = @(j) add_new2(view1,[elem(j),view1.elem_num(dim)],dim);
    replace2(view1,1,add_fun,2,dim);
end

