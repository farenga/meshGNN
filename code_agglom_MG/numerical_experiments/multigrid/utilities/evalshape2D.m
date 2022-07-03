function [dphiq,Grad]= evalshape2D(femregion, i_element, Nodes)

N = femregion.fem;
if all(femregion.nedges==3)
    loc_coord = femregion.coord(femregion.connectivity{i_element},:);
    x0=loc_coord(1,1);   % x-coordinates of vertices
    x1=loc_coord(2,1);
    x2=loc_coord(3,1);

    y0=loc_coord(1,2);   % y-coordinates of vertices
    y1=loc_coord(2,2);
    y2=loc_coord(3,2);

    BJ = [x1-x0,x2-x0;y1-y0,y2-y0];       % Jacobian of elemental map

    BJ_inv = inv(BJ);

    trans=[x0;y0];                      % translation vector
else
    BBox = femregion.BBox(i_element,:);
    x1B = BBox(1); x2B = BBox(2);
    y1B = BBox(3); y2B = BBox(4);
    
    x0=x1B;   % x-coordinates of vertices
    x1=x2B;
    x2=x2B;
    x3=x1B;
    
    y0=y1B;   % y-coordinates of vertices
    y1=y1B;
    y2=y2B;
    y3=y2B;


    %GP: this is the origin on the bouded box coordinates
    trans = (0.25) .* [ x0 + x1 + x2 + x3 ;  y0 + y1 + y2 + y3];                       % translation vector
    
    BJ = .25 .* [-x0 + x1 + x2 - x3 ,  -x0 - x1 + x2 + x3 ; -y0 + y1 + y2 - y3 , -y0 - y1 + y2 + y3];
end
% 
D = BJ(2,2)*BJ(1,1)-BJ(1,2)*BJ(2,1);

BJ_inv = [BJ(2,2) -BJ(1,2); -BJ(2,1) BJ(1,1)]/D;
trans_inv = [-BJ(2,2)*trans(1)+BJ(1,2)*trans(2);BJ(2,1)*trans(1)-BJ(1,1)*trans(2)]/D;


for k = 1:length(Nodes(:,1))
    Nodes_phys(k,:) = transpose( BJ_inv * transpose(Nodes(k,:)) + trans_inv);
end
if any(femregion.nedges==3)
    [a,b] = rstoab(Nodes_phys(:,1),Nodes_phys(:,2));
else
    a = Nodes_phys(:,1);
    b = Nodes_phys(:,2);
end

sk = 1;
for i = 0:N
  for j = 0:N-i
    if any(femregion.nedges==3)
         psi(:,sk) = Simplex2DP(a,b,i,j);
        [dpsi(:,sk,1), dpsi(:,sk,2)] = GradSimplex2DP(a,b,i,j);
        sk = sk+1;
    else
        psi(:,sk) = Hypercube2DP(a,b,i,j);
        [dpsi(:,sk,1), dpsi(:,sk,2)] = GradHypercube2DP(a,b,i,j);
        sk = sk+1;
    end
  end
end


dphiq = psi;
Grad(:,1,:)= dpsi(:,:,1);
Grad(:,2,:)= dpsi(:,:,2);
for j = 1:size(Grad,3)
    for i = 1:size(Grad,1)
        Grad(i,:,j)= Grad(i,:,j)*BJ_inv;
    end
end



