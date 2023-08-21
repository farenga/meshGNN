function [I,Iface] = IntegralOverPolygonAndFaces(Cf, i, j, v, toll)
% Evaluate Integral of the homogeneous function of degree q,
% on a polygon \Omega.
% Input data:
%  -  f homogeneous polynomial function of order q to be integrated, it is
%     f(x,y) = Cf * x^i * y^j, the degree of homogeneity is then q=i+j;
%  -  v = [x1 x2 ... xn; y1 y2 ... yn] list of polygon's verteces.
% Output:
%  -  I = integral_{\Omega} f(x,y)*dx*dy.

q = i + j; % degree of homogeneity
d = size(v,1);
m = size(v,2);

I = 0;
Iface = zeros(1,m);
for k = 1:m % loop over faces of the polygon
    
    % identify vertices of the face Fk
    if k < m
        x1 = v(1,k); x2 = v(1,k+1);
        y1 = v(2,k); y2 = v(2,k+1);
    else
        x1 = v(1,k); x2 = v(1,1);
        y1 = v(2,k); y2 = v(2,1);
    end
    
    % define ak and bk for each face Fk: they define the hyperplane Hk in which
    % the face Fk lies, that is: 
    %   the vector x=[x1 x2] \in Hk <==> <ak,x>==bk, or ak(1)*x(1) + ak(2)*x(2) == bk.
    % We choose ak such that ||ak|| = 1, that is ak will be the normal vector
    % of the face Fk.
    
    ak = (1 / sqrt( (y2-y1)^2 + (x1-x2)^2 )) .* [y2-y1 ; x1-x2];
    bk = ((y2-y1)*x1 + (x1-x2)*y1) / sqrt( (y2-y1)^2 + (x1-x2)^2 );
    
    % Evaluate the Face integral
    % smart choice of x0, but it leads to numerical errors!
    if abs(ak(1)) > toll  &&  abs(ak(2)) > toll
        if i > j
            x0(1) = 0;
            x0(2) = bk/ak(2);
        else
            x0(1) = bk/ak(1);
            x0(2) = 0;
        end
    elseif abs(ak(2)) <= toll
        x0(1) = bk/ak(1);
        x0(2) = 0;
    elseif abs(ak(1)) <= toll
        x0(1) = 0;
        x0(2) = bk/ak(2);
    end
    
    g = ReductionIntegration(Cf,i,j,x1,y1,x2,y2,x0,d);
    
    Iface(k) = g;
    I = I + bk*g;
end
I = I/(d+q);

end