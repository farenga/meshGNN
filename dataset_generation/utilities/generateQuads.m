function [mesh] = generateQuads(max_edge_q)
% grid quads
mesh = polygon(4);
for i = 1:max_edge_q
    view = RView(mesh);
    for e = eye(2)
        mesh_tmp = copy(view);
        move(mesh_tmp,e');
        merge(view,mesh_tmp);
    end
    scale(mesh,0.5);
end
end