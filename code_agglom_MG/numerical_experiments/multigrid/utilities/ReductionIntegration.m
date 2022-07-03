function g = ReductionIntegration(Cf,i,j,x1,y1,x2,y2,x0,d)
% This function is the RHS of the formula (10) described in section 3.1 of the paper
% "Numerical integration of homogeneous functions on convex and nonconvex polygons and polyhedra",
% Eric B. Chin and Jean B. Lasserre and N. Sukumar.
% 
% It is employed for computing integral_{[x1,y1]--[x2,y2]} Cf * x^i * y^j dsigma(x,y)

g = 0;
q = i + j;


% Integration of the gradient of the homogeneous function 
% that is the second term of the RHS of formula (10) in section 3.1 of the
% paper "Numerical integration of homogeneous functions on convex and nonconvex polygons and polyhedra",
% Eric B. Chin and Jean B. Lasserre and N. Sukumar.
if x0(1)*Cf*i~=0
    g = g + ReductionIntegration(x0(1)*Cf*i,i-1,j,x1,y1,x2,y2,x0,d);
end

if x0(2)*Cf*j~=0
    g = g + ReductionIntegration(x0(2)*Cf*j,i,j-1,x1,y1,x2,y2,x0,d);
end


% Now we evaluate the first term of the RHS of formula (10) in section 3.1 of the
% paper "Numerical integration of homogeneous functions on convex and nonconvex polygons and polyhedra",
% Eric B. Chin and Jean B. Lasserre and N. Sukumar.
n1 = (1/sqrt((x1-x2)^2 + (y1-y2)^2)).*[x1-x2; y1-y2];

% first intersection
dij1 = (x1-x0(1))*n1(1) + (y1-x0(2))*n1(2);
g = g + dij1*(Cf*x1^i*y1^j);

% second intersection (here n2 = -n1)
dij2 = (x0(1)-x2)*n1(1) + (x0(2)-y2)*n1(2);
g = g + dij2*(Cf*x2^i*y2^j);

g = g/(d+q-1);

end