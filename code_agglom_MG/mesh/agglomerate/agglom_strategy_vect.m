function [new_elem_list,new_weight_list] = ...
    agglom_strategy_vect(mesh,h_bar_list,aggl_fun,weights)

if mesh.elem_num(end) ==0
    error('empty mesh')
end

if nargin == 3
    weights = zeros(1,mesh.elem_num(end));
    for i = 1:mesh.elem_num(end)
        elem = RView(mesh,i,mesh.dim);
        if mesh.dim == 2
            polygon = edges2polygon(vertices(elem),edges(elem));
            P = polyshape(polygon);
            weights(i) = area(P);
        else % mesh.dim == 3
            V = vertices(elem);
            [~,weights(i)] = convhull(V);
        end
    end
end

N_h = length(h_bar_list);

if mesh.elem_num(end) <= 1
    new_elem_list(1:N_h) = {{1:mesh.elem_num(end)}};
    new_weight_list(1:N_h) = {sum(weights)};
    return
end

A = agglomerate(mesh,{1:mesh.elem_num(end)});
h = meshsize(A);
h_bar_list = sort(h_bar_list,'descend');

II = h <= h_bar_list;
new_elem_list = cell(1,N_h);
new_weight_list = cell(1,N_h);
new_elem_list(II) = {{1:mesh.elem_num(end)}};
new_weight_list(II) = {sum(weights)};

if all(II)
    return
end

if length(weights) == 2
    IDX = [1 2];
else
    IDX = aggl_fun(mesh,weights);
    IDX = fix_partition(mesh,IDX,weights);
end

% G = connectivity(mesh);
% V = vertices(G);
% plot(mesh,2)
% hold on
% plot(V(IDX==1,1),V(IDX==1,2),'o')


if all(IDX==IDX(1))
    new_elem_list(1:N_h) = {{1:mesh.elem_num(end)}};
    new_weight_list(1:N_h) = {sum(weights)};
    return
end


N_new = 2;
new_elem = cell(1,N_new);
new_weight = zeros(1,N_new);
for i = 1:N_new
    new_elem{i} = find(IDX==i);
    new_weight(i) = sum(weights(new_elem{i}));
end


M1 = RView(mesh,new_elem{1},mesh.dim);
M2 = RView(mesh,new_elem{2},mesh.dim);
W1 = weights(new_elem{1});
W2 = weights(new_elem{2});

[new_elem_list1,new_weight_list1] = agglom_strategy_vect(M1,h_bar_list(not(II)),aggl_fun,W1);
[new_elem_list2,new_weight_list2] = agglom_strategy_vect(M2,h_bar_list(not(II)),aggl_fun,W2);

N_h_left = length(new_elem_list1);
new_elem_list_left =cell(1,N_h_left);
new_weight_list_left = cell(1,N_h_left);

for i = 1:N_h_left
    for j = 1:length(new_elem_list1{i})
        new_elem_list1{i}{j} = M1.loc2glob{mesh.dim+1}(new_elem_list1{i}{j});
    end
    for j = 1:length(new_elem_list2{i})
        new_elem_list2{i}{j} = M2.loc2glob{mesh.dim+1}(new_elem_list2{i}{j});
    end
    new_elem_list_left{i} = [new_elem_list1{i},new_elem_list2{i}];
    new_weight_list_left{i} = [new_weight_list1{i},new_weight_list2{i}];
end

new_elem_list(not(II)) =  new_elem_list_left;
new_weight_list(not(II)) =  new_weight_list_left;

end