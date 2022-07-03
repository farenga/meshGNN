function [UF,CR] = quality(mesh)
dim = manifold(mesh);
diam = zeros(1,mesh.elem_num(dim+1));
for i = 1:mesh.elem_num(dim+1)
    volume = RView(mesh,i,dim);
    diam(i) = RM_incribed_cirle(volume);
end
[h,hvect] = meshsize(mesh);
CR = diam./hvect;
UF = hvect/h;
end

function diam = RM_incribed_cirle(view1)
center = conv_centr(vertices(view1),edges(view1));
if not(is_inside(view1,center))
    diam = 0;
    return
end
diam = inf;
dim = manifold(view1);
for i = 1:view1.elem_num(dim)
    dist = dist2conv(center,vertices(RView(view1,i,dim-1)));
    diam = min(diam,dist);
end
diam = 2*diam;
end

function dist = dist2conv(point,vertices)
n_vert = size(vertices,1);
C = vertices';
d = point';
Aeq = ones(1,n_vert);
beq = 1;
ub = ones(n_vert,1);
lb = zeros(n_vert,1);
options = optimoptions('lsqlin','Display','none');
x = lsqlin(C,d,[],[],Aeq,beq,lb,ub,[],options);
dist = norm(C*x-d);
end
