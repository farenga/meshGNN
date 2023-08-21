function [P] = LegendreP(x,N)

% function [P] = LegendreP(x,N)
% Purpose: Evaluate Legendre Polynomial at points x for order N and returns P[1:length(xp))]

% Turn points into row if needed.
xp = x; dims = size(xp);
if (dims(2)==1) xp = xp'; end;

PL = zeros(N+1,length(xp)); 

% Initial values P_0(x) and P_1(x)
PL(1,:) = ones(size(xp));
if (N==0) P=PL'; return; end;
PL(2,:) = xp;
if (N==1) P=PL(N+1,:)'; return; end;

for i=1:N-1
  PL(i+2,:) = 1/(i+1)*( (2*i+1)*xp.*PL(i+1,:) - i*PL(i,:));
end;

P = PL(N+1,:)';
return