clear, clc, close all
rng('default')

%% grid sorted triangles
[X,Y] = meshgrid(linspace(0,1,9));
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
figure, plot(mesh)
save([path2('grids'),'fine2D_tria'],'mesh')
%% grid random triangles
P = rand(128,2);
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
figure, plot(mesh)
save([path2('grids'),'fine2D_rand_tria'],'mesh')

%% grid random voronoi
quad = polygon(4);
mesh = voromesh2D(rand(256,2));
figure, plot(mesh)
save([path2('grids'),'fine2D_rand_voro'],'mesh')

%% grid quads
mesh = polygon(4);
for i = 1:4
    view = RView(mesh);
    for e = eye(2)
        mesh_tmp = copy(view);
        move(mesh_tmp,e');
        merge(view,mesh_tmp);
    end
    scale(mesh,0.5); % [0,2]^2 --> [0,1]^2
end
figure, plot(mesh)
save([path2('grids'),'fine2D_quads'],'mesh')

