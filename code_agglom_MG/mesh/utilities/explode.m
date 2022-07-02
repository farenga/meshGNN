function explode(mesh,numbers_dim)
hold on
c = conv_centr(vertices(mesh));
dim = manifold(mesh);
for i = 1:mesh.elem_num(dim+1)
    view = RView(mesh,i,dim);
    elem_mesh = copy(view);
    ci = conv_centr(vertices(elem_mesh));
    move(elem_mesh,2*(ci-c));    
    if nargin == 2
        plot(elem_mesh,numbers_dim)
    else
        plot(elem_mesh)
    end
end
hold off
end