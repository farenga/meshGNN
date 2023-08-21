% for each face compute the jacobian of the mapping from the refernce
% points to the pphisical point

%GP: probably this is a function for 3D code, where the hypothesis that
%faces are simplex is assumed, that is each face is a triangular. This
%function works also in 2D case, using as imput a triangular degenerated in
%the line segment representing the face

function [pphys_1D] = get_jacobian_physical_points_faces(loc_coord, node_1D)

nqn_1D=length(node_1D);
nfaces=length(loc_coord(:,1));
nodes_face =zeros(nqn_1D,2,nfaces);

x0=loc_coord(1,1);   % x-coordinates of vertices
x1=loc_coord(2,1);
x2=loc_coord(3,1);

y0=loc_coord(1,2);   % y-coordinates of vertices
y1=loc_coord(2,2);
y2=loc_coord(3,2);

nodes_face=[node_1D;zeros(1,nqn_1D)]';

BJ_face = [x1-x0,x2-x0;y1-y0,y2-y0];       % Jacobian of elemental map
trans=[x0;y0];                      % translation vector

for k=1:nqn_1D
    pphys_1D(k,:)=transpose((BJ_face*transpose(nodes_face(k,:))+trans));
end

end


