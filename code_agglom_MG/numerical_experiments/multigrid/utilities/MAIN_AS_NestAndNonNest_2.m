%==========================================
%  Main of Multigrid with ASM as Smoother
%==========================================
clear all;
close all;

n_levels = [2 3 4 5 6];
nmax = 10000;
tol = 1e-8;
smooth = [3 5];

range_p = [1];
range_h = [512 1024 2048 4096];
melem = 8; % number of elements of the coarsest grid in case nested=5 below
%range_p = [3 4];
%range_h = [512 1024 2048 4096];

% These lines come from a previous code
%filename = 'Tau_region.mat';
%load(filename,'Tau','neighbour');
%Tau_stored = Tau; clear Tau;
%neighbour_stored = neighbour; clear neighbour;

DoPureMG = 1;
DoMGPCG = 0;
DoASPCG = 0;

nested = 5;
% Legend:
% nested = 0 for non-nested polygonal Voronoi grids;
% nested = 1 for nested polygonal grids, the finest ones are Voronoi
% nested = 2 for nested polygonal grids, the finest ones are unstructured triangular
% nested = 3 for nested structured triangular grids
% nested = 4 for nested structured quadrilateral grids
% nested = 5 for nested refined grids
% nested = -1 for non-nested airfoil
% nested = -2 non-nested irregular grids

choice = 'c';
%Legend:
% a) Matrices for ASPCG and Multigrid defined as A_H = R_0^T*A_h*R_0
% b) Multigrid matrices computed while define A_H = R_0^T*A_h*R_0 only for ASPCG
% c) Multigrid matrices computed while define A_H = P_0^T*A_h*P_0 only for ASPCG,
%    where P_0 is the projection on the coarsest level
% d) Multigrid matrices inherited while define A_H = P_0^T*A_h*P_0 only for ASPCG,
%    where P_0 is the projection on the coarsest level

LastMesh = 4; % it should be 3 <= LastMesh <= 6, and it makes sense only with nested=3 or 4
mesh = 'Airfoil_h_0p25.mat';

Dati = dati_2;
A_h = []; A_H = []; P_0 = [];
path = 'MG_RESULTS/';
filename1 = ['Nest',num2str(nested)];



for choice = ['c','d']
    if choice == 'c'
        %print info
        fprintf('\n\n');
        X = '************************************************************************';
        disp(X);
        X = ['************************  MG INHERITED          ************************'];
        disp(X);
        X = '************************************************************************';
        disp(X);
        fprintf('\n');
    else
        %print info
        fprintf('\n\n');
        X = '************************************************************************';
        disp(X);
        X = ['************************ MG NON-INHERITED          ************************'];
        disp(X);
        X = '************************************************************************';
        disp(X);
        fprintf('\n');
    end
    
%     for melem = [8 16]
%         if melem == 8
%             %print info
%             fprintf('\n\n');
%             X = ['    **************  coarsest 8 elem  ************************'];
%             disp(X);
%             fprintf('\n');
%         else
%             %print info
%             fprintf('\n\n');
%             X = ['    **************  coarsest 16 elem  ************************'];
%             disp(X);
%             fprintf('\n');
%             fprintf('\n');
%         end
        
        % Generate all the grids needed for the numerical tests
        switch nested
            
            case 1
                
                % Polygonal Voronoi and agglomerated grids
                % for each fine mesh you have 4 nested grids to be loaded
                % this load is performed later
                range_h = [512 1024 2048 4096];
                %range_n_set = [1 1 1 1];
                
            case 2
                
                % Unstructured triangles and agglomerated grids
                % for each fine mesh you have 4 nested grids to be loaded
                % this load is performed later
                range_h = [582 1086 2198 4318];
                
            case 3
                
                Dati.type_mesh = 'TRIA_S';
                Dati.method = 'LDG';
                Dati.basis = 'modal';
                Taus = cell(1,LastMesh+1);
                for k = 0:LastMesh
                    melem = 2^(2*k+1);
                    X = ['    I am generating the triangular mesh with ',num2str(melem),...
                        ' elements...']; disp(X);
                    Taus{k+1} = generate_mesh_2(Dati,k);
                end
                
                Taus_all = cell(length(range_h),n_levels(end));
                index = 4;
                range_h = 2.^[7:2:2*LastMesh+1];
                for k = 1:length(range_h)
                    for l = 1:n_levels(end)
                        Taus_all{k,l} = Taus{index-l+1};
                    end
                    index = index + 1;
                end
                clear Taus;
                
            case 4
                
                Dati.type_mesh = 'QUAD';
                Dati.method = 'LDG';
                Dati.basis = 'modal';
                Taus = cell(1,LastMesh+1);
                for k = 0:LastMesh
                    melem = 4^k;
                    X = ['    I am generating the quadrilateral mesh with ',num2str(melem),...
                        ' elements...']; disp(X);
                    Taus{k+1} = generate_mesh_2(Dati,k);
                end
                
                Taus_all = cell(length(range_h),n_levels(end));
                index = 4;
                range_h = 2.^[6:2:2*LastMesh];
                for k = 1:length(range_h)
                    for l = 1:n_levels(end)
                        Taus_all{k,l} = Taus{index-l+1};
                    end
                    index = index + 1;
                end
                clear Taus;
                
            case 5
                
                % Grids defined by successive refinement
                %addpath(genpath('/Users/mac/Desktop/phd/POLYHEDRA/AGGLOMERATE'));
                %addpath(genpath('/Users/mac/Desktop/phd/POLYHEDRA/AGGLOMERATE/metis-5.0.2/metismex'));
                %melem = 8;
                X = ['    I am generating the mesh with ',num2str(melem),...
                    ' polygonal element...']; disp(X);
                Taus_all{1,1} = generate_mesh(@MbbDomain,melem,100);
                %[neighbour] = neighbours(Taus_all{1,1});
                %[error] = Plot_PolyMesh(Taus_all{1,1},1,'k',neighbour,'-','r');
                Taus_all{1,1}.tag = ones(Taus_all{1,1}.ne,1);
                nref = 4;
                for ref=1:nref
                    X = ['    I am generating the fine mesh ',num2str(ref)]; disp(X);
                    Taus_all{1,ref+1} = C_refine_grid_into_quad(Taus_all{1,ref});
                    %[neighbour] = neighbours(Taus_all{1,ref+1});
                    %[error] = Plot_PolyMesh(Taus_all{1,ref+1},ref+1,'k',neighbour,'-','r');
                end
                %[Tau_fine] = agglomerate_different_materials(Tau_fine,Tau_coarse.ne);
                %[neighbour] = neighbours(Taus_all{1,nref+1});
                %[error] = Plot_PolyMesh(Taus_all{1,nref+1},nref+1,'k',neighbour,'-','r');
                range_h = Taus_all{1,ref+1}.ne;
                for k=1:nref+1
                    Tau{k} = Taus_all{1,nref+1-k+1};
                end
                Taus_all = Tau; clear Tau;
                range_h = range_h(end:-1:1);
                a=1;
                
            case -1
                
                % Try with Tau_region_airfoil_2.mat,...,Tau_region_airfoil_6.mat
                %filename_airfoil = '/Users/mac/Desktop/ALL_PHD_MATLAB_CODES/Airfoil/distmesh_WithMetis_MOX30/Tau_region_airfoil_5.mat';
                
                % Try with Airfoil_h_0p25.mat, Airfoil_h_0p125.mat or Airfoil_h_0p0625.mat
                filename_airfoil = ['/u/pennesi/Desktop/Multilevel_NonNested_QuadratureFree/Airfoil_poly_grids/',mesh];
                load(filename_airfoil);
                range_h = [Tau_poly{1}.ne];
                
            case -2
                
                % non nested irregular grids to be loaded from the nested ones
                range_h = [2048];
                
            otherwise
                
                % Polygonal non-nested Voronoi grids
                %addpath(genpath('/u/pennesi/Desktop/PolygonClipper/'));
                Taus_all = cell(length(range_h),n_levels(end));
                for k = 1:length(range_h)
                    range = zeros(1,n_levels(end));
                    range(1) = range_h(k);
                    for l = 2:n_levels(end)
                        range(l) = range(l-1)/4;
                    end
                    
                    fprintf('\n');
                    for l = 1:n_levels(end)
                        melem = range(l);
                        X = ['    I am generating the mesh with ',num2str(melem),...
                            ' polygonal elements...']; disp(X);
                        Taus_all{k,l} = generate_mesh(@MbbDomain,melem,100);
                    end
                end
                
        end
        
        
        for p = range_p
            
            Dati = dati_2;
            Dati.fem = p;
            Dati.nqn = 2*Dati.fem + 1;
            
            [P,dP] = get_coeff_Legendre1D(p);
            
            %print info
            fprintf('\n\n');
            X = '************************************************************************';
            disp(X);
            X = ['************************ Simulations with p = ',num2str(p),' ************************'];
            disp(X);
            X = '************************************************************************';
            disp(X);
            fprintf('\n');
            
            for idx = 1:length(range_h)
                
                %print info
                fprintf('\n');
                X = ['  Cardinality of Tau_h (finest mesh): ',num2str(range_h(idx))];
                disp(X);
                
                %range = range_h(1:idx); range = range(end:-1:1);
                %EndOfLoop = idx;
                range(1) = range_h(idx);
                for k = 2:n_levels(end)
                    range(k) = range(k-1)/4;
                end
                
                switch nested
                    
                    case 1
                        
                        n_set = 1;%range_n_set(idx);
                        X = '    I am loading the set of nested polygonal grids...'; disp(X);
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
                        
                        Taus = cell(1,n_file);
                        for k = 1:n_file
                            range(k) = region_vec{k}.ne;
                            %Taus{k} = region_vec{n_file-k+1};
                            Taus{k} = region_vec{k};
                        end
                        clear region_vec; clear neigh_temp;
                        
                    case 2
                        
                        n_set = 1;
                        X = '    I am loading the set of nested unstructured triangular grids...'; disp(X);
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
                        
                        for k = 1:n_file
                            range(k) = region_vec{k}.ne;
                            %Taus{k} = region_vec{n_file-k+1};
                            Taus{k} = region_vec{k};
                        end
                        clear region_vec; clear neigh_temp;
                        
                    case -1
                        
                        for l = 1:n_levels(end)
                            Taus{l} = Tau_poly{l};
                            range(l) = Taus{l}.ne;
                        end
                        
                    case -2
                        
                        range_h = [512 1024 2048 4096];
                        n_set = 1;%range_n_set(idx);
                        X = '    I am loading the set of non-nested polygonal grids...'; disp(X);
                        n_file = 5;
                        region_vec = cell(n_file,1);
                        neigh_vec = cell(n_file,1);
                        
                        set_str = strcat('_MGG_',num2str(range_h(3)),'_set',num2str(n_set));
                        region_str = strcat('SETS/region',set_str,'.mat');
                        neigh_str = strcat('SETS/neighbour',set_str,'.mat');
                        region_temp = load(region_str);
                        temp_arg = fieldnames(region_temp);
                        region_vec{1} = region_temp.(temp_arg{1});
                        neigh_temp = load(neigh_str);
                        temp_arg = fieldnames(neigh_temp);
                        neigh_vec{1} = neigh_temp.(temp_arg{1});
                        clear region_temp neigh_temp temp_arg region_str neigh_str
                        
                        for k=4:-1:1
                            set_str = strcat('_MGG_',num2str(range_h(k)),'_set',num2str(n_set));
                            
                            i=2;
                            %for i = 0:n_file-1
                            %if i == 0
                            %    region_str = strcat('SETS/region',set_str,'.mat');
                            %    neigh_str = strcat('SETS/neighbour',set_str,'.mat');
                            %else
                            region_str = strcat('SETS/region',num2str(i),set_str,'.mat');
                            neigh_str = strcat('SETS/neighbour',num2str(i),set_str,'.mat');
                            %end
                            
                            region_temp = load(region_str);
                            temp_arg = fieldnames(region_temp);
                            region_vec{4-k+2} = region_temp.(temp_arg{1});
                            neigh_temp = load(neigh_str);
                            temp_arg = fieldnames(neigh_temp);
                            neigh_vec{4-k+2} = neigh_temp.(temp_arg{1});
                            clear region_temp neigh_temp temp_arg region_str neigh_str
                            %end
                        end
                        
                        Taus = cell(1,n_file);
                        for k = 1:n_file
                            range(k) = region_vec{k}.ne;
                            %Taus{k} = region_vec{n_file-k+1};
                            Taus{k} = region_vec{k};
                        end
                        clear region_vec; clear neigh_temp;
                        
                    otherwise
                        
                        for l = 1:n_levels(end)
                            Taus{l} = Taus_all{idx,l};
                            range(l) = Taus_all{idx,l}.ne;
                        end
                        
                end
                
                
                %%%%%%%%%%%%%%%%%%%%%
                % GENERATE MATRICES %
                %%%%%%%%%%%%%%%%%%%%%
                for k = 1:n_levels(end)
                    
                    melem = range(k);
                    
                    %print info
                    X = ['    I am generating the space with ',num2str(melem),...
                        ' elements...'];
                    disp(X);
                    
                    femregion{k} = create_dof(Dati,Taus{k});
                    neighbour = neighbours(Taus{k});
                    Matrices = QuadFreeMatrix2D_2(femregion{k},neighbour,P,dP,Dati);
                    rhs{k} = Matrices.f;
                    M{k} = Matrices.M;
                    A{k} = Matrices.A;
                    
                end
                clear Matrices, clear neighbour;
                f = rhs{1}; clear rhs;
                
                for levels = n_levels
                    
                    fprintf('\n');
                    X = ['    ==> ',num2str(levels),' LEVELS:'];
                    disp(X);
                    
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
                        
                        M_hH = Mixed_MassMatrix_QuadFree(femregion{k}, femregion{k-1}, P, p);
                        
                        %Prolongation operator R_h{k-1} : femregion{k} --> femregion{k-1}
                        R_h{k-1} = M_inv * M_hH;
                        
                        switch choice
                            
                            case 'a'
                                A{k} = (R_h{k-1}')*A{k-1}*R_h{k-1};
                                local_solvers{k-1} = MY_LOCAL_SOLVERS(femregion{k-1}, A{k-1});
                                coarse_solver{k-1} = MY_COARSE_SOLVER(femregion{k}, A{k}, R_h{k-1}, femregion{k-1}.ne);
                                
                            case 'b'
                                A_H = (R_h{k-1}')*A_h*R_h{k-1};
                                local_solvers{k-1} = MY_LOCAL_SOLVERS(femregion{k-1}, A{k-1});
                                coarse_solver{k-1} = MY_COARSE_SOLVER(femregion{k}, A_H, R_h{k-1}, femregion{k-1}.ne);
                                A_h = A_H;
                                
                            case 'c'
                                M_hH = Mixed_MassMatrix_QuadFree(femregion{levels}, femregion{k-1}, P, p);
                                P_0 = M_inv * M_hH;
                                A_H = (P_0')*A_h*P_0;
                                local_solvers{k-1} = MY_LOCAL_SOLVERS(femregion{k-1}, A{k-1});
                                coarse_solver{k-1} = MY_COARSE_SOLVER(femregion{levels}, A_H, P_0, femregion{k-1}.ne);
                                A{k} = (R_h{k-1}')*A_h*R_h{k-1};
                                A_h = A{k};
                                
                            case 'd'
                                M_hH = Mixed_MassMatrix_QuadFree(femregion{levels}, femregion{k-1}, P, p);
                                P_0 = M_inv * M_hH;
                                A_H = (P_0')*A_h*P_0;
                                local_solvers{k-1} = MY_LOCAL_SOLVERS(femregion{k-1}, A{k-1});
                                coarse_solver{k-1} = MY_COARSE_SOLVER(femregion{levels}, A_H, P_0, femregion{k-1}.ne);
                                A_h = A{k};
                        end
                    end
                    telapse = toc(tstart);
                    clear A_h; clear A_H; clear P_0;
                    
                    %print info
                    X = ['     Time taken to compute the projector: ',num2str(telapse),'.'];
                    disp(X);
                    
                    %Measure the number of ASPCG iterations
                    if(levels==n_levels(1) && DoASPCG==1)
                        z0 = zeros(femregion{1}.ndof,1);
                        tstart = tic;
                        [uh,~,iter_pcg(p,idx),~,~,~] = my_pcg(A{1}, z0, f, 100000, tol, local_solvers{1}, coarse_solver{1});
                        telapse = toc(tstart);
                        iter = iter_pcg(p,idx);
                        X = ['     - PCG completed in ',...
                            num2str(telapse),' with iter=',num2str(iter)];
                        disp(X);
                    end
                    
                    m_idx = 1;
                    for m = smooth
                        
                        fnrm2 = norm(f); r_0 = norm(f);
                        z0 = zeros(femregion{1}.ndof,1);
                        check = norm(f)/fnrm2;
                        iter = 0;
                        
                        if DoPureMG == 1
                            
                            %--------------------%
                            %     V - CYCLE      %
                            %--------------------%
                            print_info = 1;
                            tstart = tic;
                            while ( iter < nmax ) && ( check > tol )
                                if check > 100
                                    print_info = 0;
                                    disp(['     - MG with m=',num2str(m),' DOES NOT CONVERGE!']);
                                    iter = 100000;
                                else
                                    uh = W_CYCLE_AS(levels,1,A,f,z0,m,m,R_h,local_solvers,coarse_solver);
                                    check = norm( f - A{1}*uh ) / fnrm2;
                                    z0 = uh;
                                    iter = iter +1;
                                end
                            end
                            telapse = toc(tstart);
                            if print_info
                                X = ['     - MG with m=',num2str(m),' completed in ',...
                                    num2str(telapse),' with iter=',num2str(iter)];
                                disp(X);
                            end
                            r_N = norm(f-A{1}*uh);
                            rate_J(p,idx,levels-1,m_idx) = exp(1/iter*log(r_N/r_0));
                            iter_mg(p,idx,levels-1,m_idx) = iter;
                            
                        end
                        
                        if DoMGPCG == 1
                            
                            %----------------------------------------------------%
                            % Multigrid as preconditioner for Conjugate Gradient %
                            %----------------------------------------------------%
                            z0 = zeros(femregion{1}.ndof,1);
                            print_info = 1;
                            tstart = tic;
                            [uh,~,iter,~,condi,~] = MG_PCG(A, z0, f, nmax, tol, levels, m, m, R_h, local_solvers, coarse_solver);
                            telapse = toc(tstart);
                            if print_info
                                %X = ['     - MG as preconditioner with m=',num2str(m),' completed in ',...
                                %    num2str(telapse),' with iter=',num2str(iter)];
                                %disp(X);
                                X = ['     - MGPCG with m=',num2str(m),' has K(P^-1*A) = ',num2str(condi),...
                                    ' and it = ',num2str(iter)];
                                disp(X);
                            end
                            r_N = norm(f-A{1}*uh);
                            rate_J_MGPCG(p,idx,levels-1,m_idx) = exp(1/iter*log(r_N/r_0));
                            iter_MGPCG(p,idx,levels-1,m_idx) = iter;
                            
                        end
                        
                        m_idx = m_idx +1;
                    end
                end
                fprintf('\n');
            end
            clear A; clear M;
        end
        clear Taus; clear Taus_all;
        
        
%    end
end
%Save info:
file = [char(path),char(filename1)];
if DoPureMG == 1
    save(file,'iter_mg','rate_J','smooth','range_p','range_h');
end

if DoMGPCG == 1
    save(file,'iter_MGPCG','rate_J_MGPCG','smooth','range_p','range_h');
end
