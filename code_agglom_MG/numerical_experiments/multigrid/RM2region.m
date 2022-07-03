function region = RM2region(mesh)

region.ne = mesh.elem_num(3);
region.coord = vertices(mesh);

for i = 1:region.ne
    elem = RView(mesh,i,2);
    [~,vert_ids] = edges2polygon(vertices(elem),edges(elem));
    region.connectivity{i} = elem.loc2glob{1}(vert_ids);
end

region = add_mesh_statistics(region);

