clear, clc, close all
load("fine2D_brain.mat")

load("METIS_LAST.mat")

N_new = 50;
tic
IDX = aggl_metis_fun(mesh,ones(1,mesh.elem_num(end)),N_new);
aggl_time = toc;
N_new = max(IDX);
new_elem_list = cell(1,N_new);
for j = 1:N_new
    new_elem_list{j} = find(IDX==j);
end
aggl_mesh(1) = {mesh};
aggl_mesh{2} = agglomerate(mesh,new_elem_list);
save([path2('grids'),'brain_metis'],'aggl_mesh','aggl_time')
plot(aggl_mesh{2})