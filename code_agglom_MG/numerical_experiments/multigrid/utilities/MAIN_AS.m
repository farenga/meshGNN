%==========================================
%  Main of Multigrid with ASM as Smoother
%==========================================
clear all;
close all;

n_levels = [2 3 4];
nmax = 10000;
tol = 1e-8;
smooth = [3 5 8];

% range_p = [1 2];
% range_h = [512 1024 2048 4096];
range_p = [1];
range_h = [512];

filename = 'Tau_region.mat';
load(filename,'Tau','neighbour');
Tau_stored = Tau; clear Tau;
%neighbour_stored = neighbour; clear neighbour;

for p = range_p
    
    %Dati are the same for the 2 regions
    Dati = dati_2;
    Dati.fem = p;
    Dati.nqn = 2*Dati.fem + 1;
    
    [P,dP] = get_coeff_Legendre1D(p);
    
    %print info
    fprintf('\n\n');
    X = ['************************************************************************'];
    disp(X);
    X = ['************************ Simulations with p = ',num2str(p),' ************************'];
    disp(X);
    X = ['************************************************************************'];
    disp(X);
    fprintf('\n');
    
    %for levels = n_levels
    %    fprintf('\n');
    %    X = ['  ==> ',num2str(levels),' LEVELS:'];
    %    disp(X);
    
    for idx = 1:length(range_h)
        
        %print info
        fprintf('\n');
        X = ['  Cardinality of Tau_h (finest mesh): ',num2str(range_h(idx))];
        disp(X);
        
        range(1) = range_h(idx);
        for k = 2:n_levels(end)
            range(k) = range(k-1)/4;
        end
        
        %build all grids and spaces
        %for k = 1:length(range)
        ind_tau = 6 + idx;
        for k = 1:n_levels(end)
            
            melem = range(k);
            
            %print info
            X = ['    I am building the space with ',num2str(melem),...
                ' elements...'];
            disp(X);
            
            %Taus{k} = generate_mesh(@MbbDomain,melem,100);
            Taus{k} = Tau_stored{ind_tau};
            femregions{k} = create_dof(Dati,Taus{k});
            neighbour = neighbours(Taus{k});
            %Matrices = SukumarMatrix2D(femregions{k},neighbour,P,dP,Dati);
            Matrices = QuadFreeMatrix2D_2(femregions{k},neighbour,P,dP,Dati);
            
            rhs{k} = Matrices.f;
            Ms{k} = Matrices.M;
            As{k} = Matrices.A;
            
            ind_tau = ind_tau - 2;
        end
        clear Matrices, clear neighbour;
        
        f = rhs{1}; clear rhs;
        
        for levels = n_levels
            fprintf('\n');
            X = ['    ==> ',num2str(levels),' LEVELS:'];
            disp(X);
            
            for k = 1:levels
                Tau{k} = Taus{k};
                femregion{k} = femregions{k};
                M{k} = Ms{k};
                A{k} = As{k};
            end
            
            %Measure the number of CG iterations
            faiCG = 0;
            if(levels==n_levels(1) && faiCG==1)
                X = ['   I am doing the CG method...'];
                disp(X);
                tstart = tic;
                [dummy,dummy,dummy,iter_cg(p,idx)] = pcg(A{1},f,1e-8,100000);
                telapse = toc(tstart);
                iter = iter_cg(p,idx);
                X = ['     - CG completed in ',...
                    num2str(telapse),' with iter=',num2str(iter)];
                disp(X);
            end
            
            %Build projectors and matrix for sublevels
            tstart = tic;
            A_h = A{1};
            for k = 2:levels
                
                %You need the inverse of the mass matrix M_h
                M_inv = inv_SIPDG_MassMatrix(M{k-1}, ...
                    femregion{k-1}.ne, femregion{k-1}.nln);
                
                %M_hH = Mixed_MassMatrix_Sukumar(femregion{k}, femregion{k-1}, P, p);
                M_hH = Mixed_MassMatrix_QuadFree(femregion{k}, femregion{k-1}, P, p);
                
                %Prolongation operator R_h{k-1} : femregion{k} --> femregion{k-1}
                R_h{k-1} = M_inv * M_hH;
                
                %A{k} = (R_h{k-1}')*A{k-1}*R_h{k-1};
                %local_solvers{k-1} = MY_LOCAL_SOLVERS(femregion{k-1}, A{k-1});
                %coarse_solver{k-1} = MY_COARSE_SOLVER(femregion{k}, A{k}, R_h{k-1}, femregion{k-1}.ne);
            
                A_H = (R_h{k-1}')*A_h*R_h{k-1};
                local_solvers{k-1} = MY_LOCAL_SOLVERS(femregion{k-1}, A_h);
                coarse_solver{k-1} = MY_COARSE_SOLVER(femregion{k}, A_H, R_h{k-1}, femregion{k-1}.ne);
                A_h = A_H;
            end
            telapse = toc(tstart);
            clear A_h; clear A_H;
            
            %print info
            X = ['     Time taken to build the projector: ',num2str(telapse),'.'];
            disp(X);
            
            %Measure the number of PCG iterations
            faiPCG = 1;
            if(levels==n_levels(1) && faiPCG==1)
                %X = ['   I am doing the PCG method...'];
                %disp(X);
                z0 = zeros(femregion{1}.ndof,1);
                tstart = tic;
                [uh,~,iter_pcg(p,idx),~,~,~] = my_pcg(A{1}, z0, f, 100000, tol, local_solvers{1}, coarse_solver{1});
                telapse = toc(tstart);
                iter = iter_pcg(p,idx);
                X = ['     - PCG completed in ',...
                    num2str(telapse),' with iter=',num2str(iter)];
                disp(X);
            end
            
            %m_idx = 1;
            for m = smooth
                
                fnrm2 = norm(f);
                z0 = zeros(femregion{1}.ndof,1);
                check = norm(f)/fnrm2;
                iter = 0;
                
                print_info = 1;
                tstart = tic;
                if levels == 2
                    
                    %--------------------%
                    %  TWO LEVEL METHOD  %
                    %--------------------%
                    
                    while ( iter < nmax ) && ( check > tol )
                        if check > 10
                            print_info = 0;
                            disp(['     - ML with m=',num2str(m),' DOES NOT CONVERGE!']);
                            iter = 100000;
                        else
                            uh = TWO_LEVELS_AS(A,f,z0,m,m,R_h,local_solvers{1},coarse_solver{1});
                            check = norm( f - A{1}*uh ) / fnrm2;
                            z0 = uh;
                            iter = iter +1;
                        end
                    end
                    
                else
                    
                    %--------------------%
                    %     V - CYCLE      %
                    %--------------------%
                    
                    while ( iter < nmax ) && ( check > tol )
                        if check > 100
                            print_info = 0;
                            disp(['     - ML with m=',num2str(m),' DOES NOT CONVERGE!']);
                            iter = 100000;
                        else
                            uh = V_CYCLE_AS(levels,1,A,f,z0,m,m,R_h,local_solvers,coarse_solver);
                            check = norm( f - A{1}*uh ) / fnrm2;
                            z0 = uh;
                            iter = iter +1;
                        end
                    end
                    
                end
                telapse = toc(tstart);
                if print_info
                    X = ['     - ML with m=',num2str(m),' completed in ',...
                        num2str(telapse),' with iter=',num2str(iter)];
                    disp(X);
                end
                %r_N = norm(f-A{1}*uh);
                %rate_J(m_idx,p,idx) = exp(1/iter*log(r_N/r_0));
                %iter_ml(m_idx,p,idx) = iter;
                %m_idx = m_idx +1;
            end
            clear femregion; clear Tau; clear A; ...
                clear V; clear I; clear S; clear M;
            
        end
        fprintf('\n');
    end
    %end
end
%Save info:
% path = 'RESULTS_NI_JACOBI/';
% filename1 = 'errorL2';
% filename2 = 'iter_cg';
% filename3 = 'iter_ml';
% filename4 = 'rate_J';
%
% %file = [char(path),char(filename1)];
% %save(file,'errorL2');
%
% file = [char(path),char(filename2)];
% save(file,'iter_cg');
%
% file = [char(path),char(filename3)];
% save(file,'iter_ml');
%
% file = [char(path),char(filename4)];
% save(file,'rate_J');
