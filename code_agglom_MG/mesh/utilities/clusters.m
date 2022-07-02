function [clust,n_clust] = clusters(graph)

clust = zeros(1,graph.elem_num(1));
n_clust = 0;

for i = 1:graph.elem_num(1)
    if clust(i) == 0
        n_clust = n_clust + 1;
        label = n_clust;
        clust = propagate(graph,label,i,clust);
    end
end


end


function clust = propagate(graph,label,id,clust)

% label > 0
% clust(id) = 0 if no label has been assigned yet

if clust(id) > 0
    return
end

clust(id) = label;

for i = 1:graph.elem_num(2)
    edge = get1(graph,i);
    if edge(1) == id 
        clust = propagate(graph,label,edge(2),clust);
    elseif edge(2) == id 
        clust = propagate(graph,label,edge(1),clust);
    end
end


end