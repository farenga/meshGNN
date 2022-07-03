function [mesh] = generateSortedTriangles(max_nodes_st)
% grid sorted triangles
    [X,Y] = meshgrid(linspace(0,1,max_nodes_st));
    P = [X(:),Y(:)];
    DT = delaunay(P);
    mesh = RMesh(2);
    for i = 1:size(P,1)
        add_new0(mesh,P(i,:));
    end
    view = RView(mesh);
    for i = 1:length(DT)
        l1 = add1(view,DT(i,[1,2]));
        l2 = add1(view,DT(i,[2,3]));
        l3 = add1(view,DT(i,[1,3]));
        add_new2(view,[l1,l2,l3],2);
    end
end