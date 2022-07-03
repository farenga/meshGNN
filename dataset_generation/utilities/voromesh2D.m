function mesh = voromesh2D(P)
bs_ext=[0 1 1 0;0 0 1 1]';
X = P(:,1);
Y = P(:,2);
[V,C,~]=VoronoiLimit(X,Y,'bs_ext',bs_ext,'figure','off');

mesh = RMesh(2);
for i = 1:size(V,1)
    add_new0(mesh,V(i,:));
end

view = RView(mesh);
for i = 1:length(C)
    add_poly2D(view,C{i}');
end

end