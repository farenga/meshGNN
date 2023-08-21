function [P] = Hypercube2DP(a,b,i,j)

% function [P] = Hypercube2DP(a,b,i,j);
% Purpose : Evaluate 2D orthonormal polynomial
%           on hypercube at (a,b) of order (i,j).

h1 = LegendreP(a,i); h2 = LegendreP(b,j);
c = sqrt((2*i+1)*(2*j+1)/4);

P = c.*h1.*h2;

return;