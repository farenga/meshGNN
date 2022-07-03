function [u] = V_CYCLE(total_levels,lev,A,f,z0,lambda,m1,m2,femregion,P_h,R_h,M_inv)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%   MULTIGRID TWO LEVEL ALGORITHM   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(length(A) > total_levels)
    X = 'Error: number of levels is different from the number of total levels.';
    disp(X);
    exit;
end

%Iniatialize
u = z0;

if(lev == total_levels)

    u = A{lev}\f;

else
    
    %Iniatialize
    u = z0;
    
    w = 1/lambda(lev);
    
    %Pre-smoothing
    for m = 1:m1
        u = u + w*M_inv{lev}*(f - A{lev}*u);
    end
    
    %Coarse grid correction
    r = P_h{lev}*(f - A{lev}*u);
    start = zeros(femregion{lev+1}.ndof,1);
    e = V_CYCLE(total_levels,lev+1,A,r,start,lambda,m1,m2,femregion,P_h,R_h,M_inv);
    u = u + R_h{lev}*e;
    
    %Post-smoothing
    for m = 1:m2
        u = u + w*M_inv{lev}*(f - A{lev}*u);
    end
    
end