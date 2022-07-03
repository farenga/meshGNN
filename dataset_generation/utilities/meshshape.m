function mesh = meshshape(view1,nvert,nquad)

vert = farthest_first(vertices(view1),nvert);
K = convhull(vert);

mesh = RMesh(3);
view = RView(mesh);
for i = 1:nvert
    add_new0(view,vert(i,:));
end
ntria = size(K,1);
for i = 1:ntria
    P1 = vert(K(i,1),:); P2 = vert(K(i,2),:); P3 = vert(K(i,3),:);
    if 1/2*norm(cross(P2-P1,P3-P1)) < view1.tol
        error('degenarate face')
    end
    id1 = add1(view,K(i,[1 2]));
    id2 = add1(view,K(i,[2 3]));
    id3 = add1(view,K(i,[3 1]));
    add_new2(view,[id1 id2 id3],2);
end
% close all
% figure
% plot(view1)
% figure
% plot(mesh,0)
if nquad > 0
    collinear = zeros(1,mesh.elem_num(2)); % quad <-> merge 2 tria <-> remove 1 edge
    tria12 = zeros(mesh.elem_num(2),2);
    for i = 1:mesh.elem_num(2)
        tria12(i,:) = get_share(mesh,i,1);
        quad = RView(mesh,tria12(i,:),2);
        [~,~,latent] = pca(vertices(quad));
        collinear(i) = latent(3);
    end
    [~,I] = sort(collinear); % I(1:nquad) edges to remove <-> quad
    taken = [];
    k = 1;
    for i = 1:nquad % add the new quads
        while any(ismember(taken,tria12(I(k),:)))
            k = k+1;            
        end
        taken = unique([taken,tria12(I(k),:)]);        
        quad = RView(mesh,tria12(I(k),:),2);
        add_new2(mesh,setdiff(quad.loc2glob{2},I(k)),2);
        k = k+1;
    end
    ids = setdiff(1:mesh.elem_num(3),taken);
    mesh = copy(RView(mesh,ids,2)); % eliminate the old tria
end
add_new2(mesh,1:mesh.elem_num(3),3);
end
