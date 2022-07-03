function [dmodedx, dmodedy] = GradHypercube2DP(x,y,id,jd)

% function [dmodedr, dmodeds] = GradSimplex2DP(a,b,id,jd)
% Purpose: Return the derivatives of the modal basis (id,jd) on the 2D simplex at (a,b).   

Lx = LegendreP(x,id);     dLx = GradLegendreP(x, id);
Ly = LegendreP(y,jd);     dLy = GradLegendreP(y, jd);

% x-derivative
c = sqrt((2*id+1)*(2*jd+1)/4);
dmodedx = c*dLx.*Ly;

% y-derivative
dmodedy = c*Lx.*dLy;

return;
