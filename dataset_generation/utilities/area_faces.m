function A = area_faces(mesh)

A = zeros(1,mesh.elem_num(3));

for i = 1:mesh.elem_num(3)
    face = RView(mesh,i,2);
    vert = local_basis(face);
    poly = edges2polygon(vert,edges(face));
    A(i) = area(polyshape(poly));
end

end