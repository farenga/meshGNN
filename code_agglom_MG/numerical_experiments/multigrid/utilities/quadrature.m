%%--------------------------------------------------------------------
% PURPOSE:
%
% This routine computes the n Gauss-Ledendre nodes and weights on the
% reference interval (0,1)  to be used for face integrals, and
% the n^2 Gauss-Ledendre nodes and weights on the reference triangle
% (0,0), (1,0), (0,1) or on the reference square (-1,1)x(-1,1) to be used for volume integrals.
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [node_1D, w_1D, node_2D, w_2D] = quadrature(n)



[node_1D,w_1D] = gauleg(0,1,n); %GP: nodes and weights of GL in (0,1)

w_2D = [];
node_2D = [];

% n nodes and weights on the 1d interval [-1,1] GP:internal??
[x,w] = gauleg(-1,1,n);

%GP: actually the output is the nodes and weights on the triangular
%reference element
for i = 1:n
    for j = 1:n
        node_2D = [node_2D; (1+x(i))./2 , (1-x(i)).*(1+x(j))./4];
        w_2D = [w_2D, (1-x(i)).*w(i).*w(j)./8];
    end
end


%--------------------------------------------------------------------
% PURPOSE:
%
% This routine computes the n Gauss-Ledendre nodes and weights
% on a given interval (a,b)
%
% Author:
% Paola Antonietti
%--------------------------------------------------------------------

function [x,w] = gauleg(a,b,n)

m = (n+1)/2;
xm = 0.5*(b+a);
xl = 0.5*(b-a);
xx = [];

for i = 1:m
    z = cos(pi*(i-0.25)/(n+0.5)); %GP: this is a GL node on [-1,1] (??)
    while 1
        p1 = 1.0;
        p2 = 0.0;
        for j = 1:n
            p3 = p2;
            p2 = p1;
            p1 =((2.0*j-1.0)*z*p2 - (j-1.0)*p3)/j;
        end
        pp = n*(z*p1-p2)/(z*z-1.0);
        z1 = z;
        z = z1-p1/pp;
        if (abs(z-z1)<eps), break, end
    end
    xx(i) = xm-xl*z;
    xx(n+1-i) = xm+xl*z;
    ww(i) = 2.0*xl/((1.0-z*z)*pp*pp);
    ww(n+1-i) = ww(i);
end

x = xx;
w = ww;