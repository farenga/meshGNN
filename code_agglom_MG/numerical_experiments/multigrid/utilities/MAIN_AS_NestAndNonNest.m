%==========================================
%  Main of Multigrid with ASM as Smoother
%==========================================
clear all;
close all;

n_levels = [2 3 4];
nmax = 10000;
tol = 1e-6;
smooth = [3 5 8];

range_p = [1 2 3 4];
range_h = [16 64 256 1024 4096]; start = 2;

filename = 'Tau_region.mat';
load(filename,'Tau','neighbour');
Tau_stored = Tau; clear Tau;
%neighbour_stored = neighbour; clear neighbour;

nested = -2;

% Legend:
% nested = 0 for non-nested polygonal Voronoi grids;
% nested = 1 for nested polygonal grids, the finest ones are Voronoi
% nested = 2 for nested structured quadrilateral grids
% nested = 3 for nested structured triangular grids
% nested = 4 for nested polygonal grids, the finest ones are unstructured triangular
% nested = -1 for non-nested unstructured triangular grids
% nested = -2 for non-nested grids, the finest ones are Voronoi polygonal, 
%     while the coarser ones are unstructured triangular

% Generate all the grids needed for the numerical tests
if nested == 1
    
    range_h = [512 1024 2048 4096]; start = 1;
    range_n_set = [1 1 1 1];
    
elseif nested == 2
    
    Dati = dati_2;
    Dati.type_mesh = 'QUAD'; Dati.method = 'LDG'; Dati.basis = 'modal';
    range_h = [4 16 64 256 1024 4096]; start = n_levels(end);
    for k = 1:length(range_h)
        melem = 4^k;
        X = ['    I am building the quadrilateral mesh with ',num2str(melem),...
            ' elements...']; disp(X);
        [Taus{k}] = generate_mesh_2(Dati,k);
    end
    
elseif nested == 3
    
    Dati = dati_2;
    Dati.type_mesh = 'TRIA_S'; Dati.method = 'LDG'; Dati.basis = 'modal';
    range_h = [8 16 64 256 1024 4096]; start = n_levels(end);
    range_h = 2*range_h;
    for k = 1:length(range_h)
        melem = 2^(2*k+1);
        X = ['    I am generating the triangular mesh with ',num2str(melem),...
            ' elements...']; disp(X);
        [Taus{k}] = generate_mesh_2(Dati,k);
    end
    
elseif nested == 4
    
    % unstructured triangles
    range_h = [582 1086 2198 4318]; start = 1;
    
    % for each fine mesh you have 4 nested grids to be loaded
    
elseif nested == -1
    
    Dati = dati_2;
    
    hmax = 0.8;
    for k = 1:length(range_h)
        X = ['    I am generating the triangular mesh with hmax=',num2str(hmax/2),...
            '...']; disp(X);
        [Taus{k}] = generate_UnstructTriMeshOnSquare(Dati,hmax);
        help(k) = Taus{k}.ne;
        hmax = hmax/2;
    end
    range_h = help;
    clear help;
    
elseif nested == -2
    
    range_h = [512 1024 2048 4096]; start = 1;
    X = ['  The finest one is polygonal, the coeaser ones are triangular.']; disp(X);
    
elseif nested == -3
    %NOT IMPLEMENTED YET
    X = ['  The finest one is polygonal, the coarser ones are structured nested quadrilateral.']; disp(X);

elseif nested == -4
    %NOT IMPLEMENTED YET
    X = ['  The finest one is polygonal, the coarser ones are structured nested triangular.']; disp(X);
    
else
    
    %addpath(genpath('/u/pennesi/Desktop/PolygonClipper/'));
    range_h = [8 16 32 64 128 256 512 1024 2048 4096]; start = 7;
    for k = 1:length(range_h)
        
        melem = range_h(k);
        X = ['    I am generating the mesh with ',num2str(melem),...
            ' ...']; disp(X);
        Taus{k} = generate_mesh(@MbbDomain,melem,100);
        
    end
    
end


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
    
    for idx = start:length(range_h)
        
        %print info
        fprintf('\n');
        X = ['  Cardinality of Tau_h (finest mesh): ',num2str(range_h(idx))];
        disp(X);
        
        %range(1) = range_h(idx);
        %for k = 2:n_levels(end)
        %    range(k) = range(k-1)/4;
        %end
        
        
        if nested == 1
            
            n_set = range_n_set(idx);
            X = ['    I am loading the set of nested polygonal grids...']; disp(X);
            n_file = 5;
            set_str = strcat('_MGG_',num2str(range_h(idx)),'_set',num2str(n_set));
            
            region_vec = cell(n_file,1);
            neigh_vec = cell(n_file,1);
            for i = 0:n_file-1
                if i == 0
                    region_str = strcat('SETS/region',set_str,'.mat');
                    neigh_str = strcat('SETS/neighbour',set_str,'.mat');
                else
                    region_str = strcat('SETS/region',num2str(i),set_str,'.mat');
                    neigh_str = strcat('SETS/neighbour',num2str(i),set_str,'.mat');
                end
                
                region_temp = load(region_str);
                temp_arg = fieldnames(region_temp);
                region_vec{i+1} = region_temp.(temp_arg{1});
                neigh_temp = load(neigh_str);
                temp_arg = fieldnames(neigh_temp);
                neigh_vec{i+1} = neigh_temp.(temp_arg{1});
                clear region_temp neigh_temp temp_arg region_str neigh_str
            end
            
            for k = 1:n_levels(end)
                range(k) = region_vec{k}.ne;
                Tau_here{k} = region_vec{k};
            end
            clear region_vec; clear neigh_temp;
            
        elseif nested == 4
            
            n_set = 1;
            X = ['    I am loading the set of nested unstructured triangular grids...']; disp(X);
            n_file = 5;
            set_str = strcat('_MGG_tria_',num2str(range_h(idx)),'_set',num2str(n_set));
            
            region_vec = cell(n_file,1);
            neigh_vec = cell(n_file,1);
            for i = 0:n_file-1
                if i == 0
                    region_str = strcat('SETS/region',set_str,'.mat');
                    neigh_str = strcat('SETS/neighbour',set_str,'.mat');
                else
                    region_str = strcat('SETS/region',num2str(i),set_str,'.mat');
                    neigh_str = strcat('SETS/neighbour',num2str(i),set_str,'.mat');
                end
                
                region_temp = load(region_str);
                temp_arg = fieldnames(region_temp);
                region_vec{i+1} = region_temp.(temp_arg{1});
                neigh_temp = load(neigh_str);
                temp_arg = fieldnames(neigh_temp);
                neigh_vec{i+1} = neigh_temp.(temp_arg{1});
                clear region_temp neigh_temp temp_arg region_str neigh_str
            end
            
            for k = 1:n_levels(end)
                range(k) = region_vec{k}.ne;
                Tau_here{k} = region_vec{k};
                %Taus{k} = region_vec{k};
            end
            clear region_vec; clear neigh_temp;
            
            EndOfLoop = 5;
            
        elseif nested == -2
            
            range = zeros(1,n_levels(end));
            range(1) = range_h(idx);
            X = ['    I am generating the mesh with ',num2str(range(1)),...
                ' elements ...']; disp(X);
            Tau_here{1} = generate_mesh(@MbbDomain,range(1),100);
            hmax = Tau_here{1}.h*2;
            for k = 2:n_levels(end)
                X = ['    I am generating the unstructured triangular mesh with hmax=',num2str(hmax),...
                    '...']; disp(X);
                [Tau_here{k}] = generate_UnstructTriMeshOnSquare(Dati,hmax);
                range(k) = Tau_here{k}.ne;
                hmax = hmax*2;
            end
            EndOfLoop = idx;
           
        elseif nested == 2 || nested == 3
            
            range = range_h(idx-n_levels(end)+1:idx); range = range(end:-1:1);
            index = 1;
            for k = [idx:-1:idx-n_levels(end)+1]
                Tau_here{index} = Taus{k};
                index = index + 1;
            end
            
        elseif nested == 0
            
            range = range_h(idx-n_levels(end):idx); range = range(end:-1:1);
            index = 1;
            for k = [idx:-2:idx-2*(n_levels(end)-1)]
                Tau_here{index} = Taus{k};
                index = index + 1;
            end
            
        end
        
        
        
        %for k = 1:length(range)
        %ind_tau = 6 + idx;
        for k = 1:n_levels(end)
            
            melem = range(k);
            
            %print info
            X = ['    I am generating the space with ',num2str(melem),...
                ' elements...'];
            disp(X);
            
            %Taus{k} = generate_mesh(@MbbDomain,melem,100);
            %Taus{k} = Tau_stored{ind_tau};
            femregions{k} = create_dof(Dati,Tau_here{k});
            neighbour = neighbours(Tau_here{k});
            Matrices = QuadFreeMatrix2D(femregions{k},neighbour,P,dP,Dati);
            
            rhs{k} = Matrices.f;
            Ms{k} = Matrices.M;
            As{k} = Matrices.A;
            
            %ind_tau = ind_tau - 2;
        end
        clear Matrices, clear neighbour;
        
        f = rhs{1}; clear rhs;
        
        for levels = n_levels
            
            fprintf('\n');
            X = ['    ==> ',num2str(levels),' LEVELS:'];
            disp(X);
            
            for k = 1:levels
                %Tau{k} = Tau_here{k};
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
            for k = 2:levels
                
                %You need the inverse of the mass matrix M_h
                M_inv = inv_SIPDG_MassMatrix(M{k-1}, ...
                    femregion{k-1}.ne, femregion{k-1}.nln);
                
                M_hH = Mixed_MassMatrix_QuadFree(femregion{k}, femregion{k-1}, P, p);
                
                %Prolongation operator R_h{k-1} : femregion{k} --> femregion{k-1}
                R_h{k-1} = M_inv * M_hH;
                
                local_solvers{k-1} = MY_LOCAL_SOLVERS(femregion{k-1}, A{k-1});
                coarse_solver{k-1} = MY_COARSE_SOLVER(femregion{k}, A{k}, R_h{k-1}, femregion{k-1}.ne);
            end
            telapse = toc(tstart);
            
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
