function coarse_mesh = agglomerate(fine_mesh,new_elem)
d = fine_mesh.dim;
coarse_mesh = RMesh(d);
view = RView(coarse_mesh);
for i = 1:length(new_elem)
	agglom = copy(boundary(RView(fine_mesh,new_elem{i},d)));
    add_new2(agglom,1:agglom.elem_num(d),d);
    merge(view,agglom);
end
end