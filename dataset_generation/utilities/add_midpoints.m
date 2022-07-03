function add_midpoints(mesh)
for i = 1:mesh.elem_num(2)
    edge = get1(mesh,i);
    c = mean(get0(mesh,edge));
    ic = add_new0(mesh,c);    
    add_fun = @(j)  add_new1(mesh,[edge(j),ic]);
    replace1(mesh,i,add_fun,2);
end
end