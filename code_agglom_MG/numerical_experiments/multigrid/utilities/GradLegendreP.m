function [dP] = GradLegendreP(x,N)

% function [dP] = GradLegendreP(x,N)
% Purpose: Evaluate Grad Legendre Polynomial at points x for order N and returns P[1:length(xp))]


% Turn points into row if needed.
xp = x; dims = size(xp);
if (dims(2)==1) xp = xp'; end;

dPL = zeros(N+1,length(xp)); 

% Initial values P_0(x) and P_1(x)
dPL(1,:) = 0.0;
if (N==0) dP=dPL'; return; end;
dPL(2,:) = 1.0;
if (N==1) dP=dPL(N+1,:)'; return; end;

vec = [N-1:-2:0];

dP = 0;
for i=vec
    [P] = LegendreP(x,i);
    dP = dP + 2*P/(2/(2*i+1));
end

return