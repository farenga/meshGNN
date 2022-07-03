function [BJ, BJ_inv, pphys_2D] = map_local_physical_points(loc_coord, nodes_2D)

%GP: map the nodes_2D of the reference triangular element ((0,0) (1,0)
%(0,1)) on the triangle ((x0,y0) (x1,y1) (x2,y2))

x0=loc_coord(1,1);   % x-coordinates of vertices
x1=loc_coord(2,1);
x2=loc_coord(3,1);

y0=loc_coord(1,2);   % y-coordinates of vertices
y1=loc_coord(2,2);
y2=loc_coord(3,2);

BJ = [x1-x0,x2-x0;y1-y0,y2-y0];       % Jacobian of elemental map

BJ_inv = inv(BJ);

trans=[x0;y0];                      % translation vector


for k=1:length(nodes_2D(:,1))
    pphys_2D(k,:) = transpose( BJ * transpose(nodes_2D(k,:)) + trans);
end
