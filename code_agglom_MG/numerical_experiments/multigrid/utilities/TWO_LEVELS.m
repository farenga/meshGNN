function [u] = TWO_LEVELS(A,f,z0,lambda,m1,m2,femregion,P_h,R_h,M_inv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   MULTIGRID TWO LEVEL ALGORITHM   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(A)>2)
    X = 'Error: number of levels is different from 2.';
    disp(X);
    exit;
end

%Iniatialize
u = z0;

lev = 1;

w = 1/lambda(1);

%Pre-smoothing: m1 steps of Richardson iteration
for m = 1:m1
    u = u + w * M_inv{1} * (f - A{1}*u);
end

%Coarse grid correction
r = P_h{1}*(f - A{1}*u);
e = A{2}\r;
u = u + R_h{1}*e;

%Post-smoothing: m2 steps of Richardson iteration
for m = 1:m2
    u = u + w * M_inv{1} * (f - A{1}*u);
end



