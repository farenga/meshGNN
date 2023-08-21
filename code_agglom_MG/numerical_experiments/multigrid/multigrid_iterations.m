function [iter,iter_cg] = multigrid_iterations(aggl_mesh,m,p)

%==========================================
%  Main of Multigrid with ASM as Smoother
%==========================================

%% parameters
nmax = 10000;
tol = 1e-6;
do_CG = true;

%% grid 
levels = length(aggl_mesh);
Tau = cell(1,levels);
Nelem = zeros(1,levels);
for i = 1:levels
    Tau{i} = RM2region(aggl_mesh{i});
    Nelem(i) = Tau{i}.ne;
end

%% Dati are the same for the 2 regions
Dati = dati_2;
Dati.fem = p;
Dati.nqn = 2*Dati.fem + 1;
% m = 3*p^2;

%% build all grids and spaces
[P,dP] = get_coeff_Legendre1D(p);
femregion = cell(1,levels);
M = cell(1,levels);
A = cell(1,levels);
for k = 1:levels
    disp(['Building the space with ',num2str(Nelem(k)),' elements...']);
    femregion{k} = create_dof(Dati,Tau{k});
    neighbour = neighbours(Tau{k});
    Matrices = QuadFreeMatrix2D_2(femregion{k},neighbour,P,dP,Dati);
    if k == 1
        f = Matrices.f;
    end
    M{k} = Matrices.M;
    A{k} = Matrices.A;   
end

%% Measure the number of CG iterations
if do_CG
    disp('I am doing the CG method...');
    tstart = tic;
    [~,~,~,iter_cg] = pcg(A{1},f,tol,nmax); % aggiungo chiamata a quello che sta sotto
    telapse = toc(tstart);
    disp(['- CG time: ',num2str(telapse),' CG iter: ',num2str(iter_cg)]);
end

%% Build projectors and matrix for sublevels
tstart = tic;
A_h = A{1};
R_h = cell(1,levels-1);
local_solvers = cell(1,levels-1);
coarse_solver = cell(1,levels-1);
for k = 2:levels
    
    % You need the inverse of the mass matrix M_h
    M_inv = inv_SIPDG_MassMatrix(M{k-1}, ...
        femregion{k-1}.ne, femregion{k-1}.nln);

    M_hH = Mixed_MassMatrix_QuadFree(femregion{k}, femregion{k-1}, P, p);
    
    % Prolongation operator R_h{k-1} : femregion{k} --> femregion{k-1}
    R_h{k-1} = M_inv * M_hH;
    
    A_H = (R_h{k-1}')*A_h*R_h{k-1};
    local_solvers{k-1} = MY_LOCAL_SOLVERS(femregion{k-1}, A_h);
    coarse_solver{k-1} = MY_COARSE_SOLVER(femregion{k}, A_H, R_h{k-1},...
        femregion{k-1}.ne);
    A_h = A_H;
end
telapse = toc(tstart);
disp(['Time taken to build the projector: ',num2str(telapse)]);

%% smoothing
norm_f = norm(f);
z0 = zeros(femregion{1}.ndof,1);
check = 1;
iter = 0;

print_info = true;
tstart = tic;

while iter < nmax && check > tol
    if check > 10
        print_info = false;
        disp(['- ML with m = ',num2str(m),' DOES NOT CONVERGE!']);
        iter = 100000;
    else
        if levels == 2
            % TWO LEVEL METHOD 
            uh = TWO_LEVELS_AS(A,f,z0,m,m,R_h,local_solvers{1},coarse_solver{1});
        else
            % V-CYCLE
            uh = V_CYCLE_AS(levels,1,A,f,z0,m,m,R_h,local_solvers,coarse_solver);
        end
        check = norm( f - A{1}*uh ) / norm_f;
        z0 = uh;
        iter = iter +1;
    end
end
telapse = toc(tstart);
if print_info
    disp(['- ML with m = ',num2str(m),' completed in ',num2str(telapse),...
        ' with iter: ',num2str(iter)]);
end

% r_N = norm(f-A{1}*uh);
% rate_J = exp(1/iter*log(r_N/r_0));
% iter_ml = iter;

end