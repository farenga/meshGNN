function [mesh] = generateRandomTriangles(rt)
% grid random triangles    
    P = rand(rt,2);
    quad = polygon(4);
    add_midpoints(quad)
    add_midpoints(quad)
    add_midpoints(quad)
    P = [P;vertices(quad)];
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