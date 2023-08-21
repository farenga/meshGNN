%===========================================================%
%  Main of Multigrid with Richardson iteration as Smoother  %
%===========================================================%

% Model problem:
%   - \delta u = f    on \Omega = (0,1)^2
%   u = 0     on \partial \Omega
%
% let f to be the forcing term such that uex = (x.^2-x).*(y.^2-y).*exp(x)
% is the exact solution of the problem


% We test the Multilevel method with non nested subspaces: here the
% subspaces are obtained by defining the Voronoi meshes independently, such
% that if N is the number of element of the finer mesh, then the coarser
% has N/4 polygonal elements.

% In order to run the code you need Polygon Clipper, see
% https://uk.mathworks.com/matlabcentral/fileexchange/8818-polygon-clipper?requestedDomain=true

clear all;
close all;


nmax = 10000;
tol = 1e-8;

% test the multilevel for different number of levels
n_levels = [3];

% test the multilevel for different number of Richardson smoothing steps
smooth = [3 5 7 9];

% test for different polynomial degrees
range_p = [1];

% number of finest mesh
range_h = [80 160 320 640];

% Use already builded grids
%filename = 'Tau_region.mat';
%load(filename,'Tau','neighbour');
%Tau_stored = Tau; clear Tau;
%neighbour_stored = neighbour; clear neighbour;

for p = range_p
    
    %Dati are the same for the 2 regions
    Dati = dati;
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
        %ind_tau = 6 + idx;
        for k = 1:n_levels(end)
            
            melem = range(k);
            
            %print info
            X = ['    I am building the space with ',num2str(melem),...
                ' elements...'];
            disp(X);
            
            Taus{k} = generate_mesh(@MbbDomain,melem,100);
            %Taus{k} = Tau_stored{ind_tau};
            femregions{k} = create_dof(Dati,Taus{k});
            neighbour = neighbours(Taus{k});
            Matrices = QuadFreeMatrix2D(femregions{k},neighbour,P,dP,Dati);
            
            rhs{k} = Matrices.f;
            Ms{k} = Matrices.M;
            As{k} = Matrices.A;
            
            %ind_tau = ind_tau - 2;
        end
        clear Matrices, clear neighbour;
        
        f = rhs{1}; clear rhs;
        
        for total_levels = n_levels
            
            fprintf('\n');
            X = ['    ==> ',num2str(total_levels),' LEVELS:'];
            disp(X);
            
            for k = 1:total_levels
                Tau{k} = Taus{k};
                femregion{k} = femregions{k};
                M{k} = Ms{k};
                A{k} = As{k};
            end
            
            %Build L2-projectors and matrix for sublevels
            tstart = tic;
            for k = 2:total_levels
                
                %You need the inverse of the mass matrix M_H
                M_inv{k-1} = inv_SIPDG_MassMatrix(M{k-1}, ...
                    femregion{k-1}.ne, femregion{k-1}.nln);
                
                M_hH = Mixed_MassMatrix_QuadFree(femregion{k}, femregion{k-1}, P, p);
                
                %Prolongation operator R_h{k-1} : femregion{k} --> femregion{k-1}
                R_h{k-1} = M_inv{k-1} * M_hH;
                P_h{k-1} = R_h{k-1}';
                
            end
            telapse = toc(tstart);
            
            %print info
            X = ['     Time taken to build the projector: ',num2str(telapse),'.'];
            disp(X);
            
            %Uper bound for Richardson smoother
            opts.disp = 0;
            opts.tol = 1e-8;
            for k = 1:total_levels-1
                lambda(k) = eigs(M_inv{k}*A{k},1,'lm',opts);
            end
            
            % loop over different smoothing steps
            for m = smooth
                
                fnrm2 = norm(f);
                z0 = zeros(femregion{1}.ndof,1);
                check = norm(f)/fnrm2;
                iter = 0;
                
                print_info = 1;
                tstart = tic;
                if total_levels == 2
                    
                    %--------------------%
                    %  TWO LEVEL METHOD  %
                    %--------------------%
                    
                    while ( iter < nmax ) && ( check > tol )
                        if check > 10
                            print_info = 0;
                            disp(['     - ML with m=',num2str(m),' DOES NOT CONVERGE!']);
                            iter = 100000;
                        else
                            uh = TWO_LEVELS(A,f,z0,lambda,m,m,femregion,P_h,R_h,M_inv);
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
                            uh = V_CYCLE(total_levels,1,A,f,z0,lambda,m,m,femregion,P_h,R_h,M_inv);
                            check = norm( f - A{1}*uh ) / fnrm2;
                            z0 = uh;
                            iter = iter +1;
                        end
                    end
                    
                end
                telapse = toc(tstart);
                if print_info
                    X = ['     - ML with m=',num2str(m),' coverges with iter=',num2str(iter)];
                    disp(X);
                end
                
            end
            clear femregion; clear Tau; clear A; ...
                clear V; clear I; clear S; clear M;
            
        end
        fprintf('\n');
    end
    %end
end

