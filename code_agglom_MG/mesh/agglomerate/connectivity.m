function [graph,W] = connectivity(mesh,dim)
% graph nodes have the same ordering of mesh elements

if nargin < 2
    dim = mesh.dim;
end

graph = RMesh(mesh.dim);

for i = 1:mesh.elem_num(dim+1)
    view = RView(mesh,i,dim);
    c = mean(vertices(view));
    add_new0(graph,c);
end

W = zeros(mesh.elem_num(dim+1),mesh.elem_num(dim+1));
if dim == 2
    area = length_edges(mesh);
else
    area = area_faces(mesh); % to write function that calculates area
end

Vgraph = RView(graph);

if dim > 0
    for i = 1:mesh.elem_num(dim) % for every sub-element
        % elements sharing that sub-elem are connected 
        connected = get_share(mesh,i,dim-1);
        n = length(connected);
        for j = 1:(n-1)
            for k = (j+1):n
                link = [connected(j),connected(k)];
                W(connected(j),connected(k)) = ...
                    W(connected(j),connected(k)) + area(i);
                W(connected(k),connected(j)) = W(connected(j),connected(k));
                add1(Vgraph,link);
            end
        end
    end 
end

end