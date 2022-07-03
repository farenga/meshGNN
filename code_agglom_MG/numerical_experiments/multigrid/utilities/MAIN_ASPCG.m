%==========================================
%  Main of AS preconditioner for CG method
%==========================================
clear all;
close all;

%n_levels = [2];
nmax = 10000;
tol = 1e-8;

% add the path to PolygonClipper
addpath(genpath('C:\Users\giorgio.pennesi\Desktop\PhD\PolygonClipper'));
        
range_p = [2 3 4 5 6 7 8];
%range_h = [16 64 256 1024]; start = 4;
range_h = [8 262]; start = 2;
range_rho = 10.^[0];

%startLine615 = 2;
startLine615 = 2;

%range_p = [1:10];
%range_h = [1024 4096]; start = 2;

%filename = 'Tau_region.mat';
%load(filename,'Tau','neighbour');
%Tau_stored = Tau; clear Tau;
%neighbour_stored = neighbour; clear neighbour;

choice_of_A = 3;
% Legend:
% choice_of_A = 1 for matrices computed on each level
% choice_of_A = 2 for A{1} computed and A{k} such that
%                  V{k} = (R_h')*V{1}*R_h;
%                  S{k} = (h/H)*(R_h')*S{1}*R_h;
%                  I{k} = (R_h')*I{1}*R_h;
%                  A{k} = V{k} - transpose(I{k}) - I{k} + S{k};
% choice_of_A = 3 for A{1} computed and A{k} = (R_h')*A{1}*R_h;

choice_of_size = 2;
% Legend:
% choice_of_size = 1 for h < H < HH, here in particular HH = 2*H;
% choice_of_size = 2 for h = H < HH;
% choice_of_size = 3 for h < H = HH;

nested = -7;
% Legend:
% nested = -7 for some particular non-nested grids, where the coarse is 
%             Voronoi, employed for the p-dependence test
% nested = -6 for some particular non-nested grids employed for
%             the p-dependence test
% nested = -5 for the airfoil domain
% nested = -4 for non-nested grids, the finest ones are Voronoi polygonal,
%             while the coarser ones are unstructured triangular
% nested = -3 for non-nested grids, the finest ones are Voronoi polygonal,
%             while the coarser ones are structured quadrilateral
% nested = -2 for non-nested grids, the finest ones are Voronoi polygonal,
%             while the coarser ones are structured triangular
% nested = -1 for non-nested unstructured triangular grids
% nested = 0 for non-nested polygonal Voronoi grids;
% nested = 1 for nested polygonal grids, the finest ones are Voronoi
% nested = 2 for nested polygonal grids, the finest ones are unstructured triangular
% nested = 3 for nested structured triangular grids
% nested = 4 for nested structured quadrilateral grids
% nested = 5 the fine is triangular, the coarse is quadrilateral, it is
%            like the Dr.-Krz. 2016 paper
% nested = 6 experiment for testing rho with Chessboard pattern on nested
%            polygons
% nested = 7 attempt with L-shaped domain, variable rho, coarse grid
%            elements are convex shaped

q = 1;
% Legend:
% q = 0 actually means that q is chosen as q = p
% q = K with k >= 1 means that q is constant 

%mesh = 'Airfoil_h_0p25.mat'; % if nested=-5 then select between h_0p25, h_0p125 and h_0p0625

Dati = dati_2;

% Generate all the grids needed for the numerical tests
switch nested
    
    case 1
        
        % Polygonal Voronoi and agglomerated grids
        
        range_h = [512 1024 2048 4096]; start = 1;
        % for each fine mesh you also have 4 nested sub-grids to be loaded
        % this load is performed later
        
    case 2
        
        % Unstructured triangles and agglomerated grids
        
        range_h = [582 1086 2198 4318]; start = 1;
        % for each fine mesh you also have 4 nested sub-grids to be loaded
        % this load is performed later
        
    case 3
        
        Dati = dati_2;
        Dati.type_mesh = 'TRIA_S';
        Dati.method = 'LDG';
        Dati.basis = 'modal';
        range_h = 2*range_h;
        for k = 2:length(range_h)+1
            
            melem = 2^(2*k+1);
            X = ['    I am generating the triangular mesh with ',num2str(melem),...
                ' elements...']; disp(X);
            [Taus{k-1}] = generate_mesh_2(Dati,k);
        end
        
    case 4
        
        Dati = dati_2;
        Dati.type_mesh = 'QUAD';
        Dati.method = 'LDG';
        Dati.basis = 'modal';
        
        index = -2 + log2(range_h(1))/2;
        for k = 2:length(range_h)+1
            melem = 4^(k+index);
            X = ['    I am generating the quadrilateral mesh with ',num2str(melem),...
                ' elements...']; disp(X);
            [Taus{k-1}] = generate_mesh_2(Dati,k+index);
        end
        
    case 5
        
        Dati = dati_2;
        Dati.type_mesh = 'QUAD';
        Dati.method = 'LDG';
        Dati.basis = 'modal';
        
        index = -2 + log2(range_h(1))/2;
        %if range_h(1) == 4
        %    index = -1;
        %elseif range_h(1) == 16
        %    index = 0;
        %elseif range_h(1) == 64
        %    index = 1;
        %elseif range_h(1) == 256
        %    index = 2;
        %elseif range_h(1) == 1024
        %    index = 3;
        %elseif range_h(1) == 4096
        %    index = 4;
        %end
        
        %for k = 2:length(range_h)+1
        k = 2+index;
        %melem = 4^(k+index);
        melem = 4^(k);
        X = ['    I am generating the quadrilateral coarse grid with ',num2str(melem),...
            ' elements...']; disp(X);
        [Taus{1}] = generate_mesh_2(Dati,k);
        %    [Taus{k-1}] = generate_mesh_2(Dati,k+index);
        %end
        
        
        Dati = dati_2;
        Dati.type_mesh = 'TRIA_S';
        Dati.method = 'LDG';
        Dati.basis = 'modal';
        %range_h = 2*range_h;
        %for k = 2:length(range_h)+1
        k = 5;
        melem = 2^(2*k+1);
        X = ['    I am generating the triangular fine grid with ',num2str(melem),...
            ' elements...']; disp(X);
        [Taus{2}] = generate_mesh_2(Dati,k);
        %[Taus{k-1}] = generate_mesh_2(Dati,k);
        %end
        
    case 6
        
        % Employed only to show the dependence on polygonal pattern of rho
        start = 2;
        idx = 1;
        range_h = [1024];
        n_set = 1;
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
        Dati.permeability = 'PolygonalChessboard';
        Dati.DiffusionCoefficientRegion = region_vec{5}; % Chessboard
        Taus{1} = region_vec{5}; range_h(1) = Taus{1}.ne;% Tau_H
        Taus{2} = region_vec{3}; range_h(2) = Taus{2}.ne;% Tau_h
        
        
        
        
        
    case 7
        
        % Polygonal nested Voronoi grids on L-shaped domain
        
        
        
        %addpath(genpath('/u/pennesi/Desktop/PolygonClipper/'));
        %addpath(genpath('/Users/mac/Desktop/phd/POLYHEDRA/AGGLOMERATE'));
        %addpath(genpath('/Users/mac/Desktop/phd/POLYHEDRA/AGGLOMERATE/metis-5.0.2/metismex'));
        addpath(genpath('C:\Users\giorgio.pennesi\Desktop\PhD\Multilevel_NonNested_QuadratureFree_MOX30\MeshIlario'));
%         melem = 16;
%         X = ['    I am generating the mesh with ',num2str(melem),...
%             ' polygonal elements on non convex doamin...']; disp(X);
%         Tau_coarse = generate_mesh_L_shaped(@MbbDomain_L_shaped,melem,100);
%         %Tau_coarse = generate_mesh(@MbbDomain_Non_Convex,melem,100);
%         [neighbour] = neighbours(Tau_coarse);
%         [error] = Plot_PolyMesh(Tau_coarse,1,'k',neighbour,'-','r');
%         
%         X = ['    I am generating the fine mesh on the same domain...']; disp(X);
%         %Tau_fine] = generate_fine_mesh_CartesianIntersection(Tau_coarse,3);
%         %[neighbour] = neighbours(Tau_fine);
%         %[error] = Plot_PolyMesh(Tau_fine,2,'k',neighbour,'-','r');
%         nref = 4; Tau_fine = Tau_coarse;
%         for ref=1:nref
%             Tau_fine = C_refine_grid_into_quad(Tau_fine);
%             %[neighbour] = neighbours(Tau_fine_here);
%             %[error] = Plot_PolyMesh(Tau_fine_here,ref+1,'k',neighbour,'-','r');
%             %Tau_temp = Tau_fine;
%         end
%         [Tau_fine] = agglomerate_different_materials(Tau_fine,Tau_coarse.ne);
%         [neighbour] = neighbours(Tau_fine);
%         [error] = Plot_PolyMesh(Tau_fine,nref+1,'k',neighbour,'-','r');
%         a=1;
        
        
        
        % Cartesian grids on L-shaped domain
        %         Dati = dati_2;
        %         Dati.type_mesh = 'QUAD';
        %         Dati.method = 'LDG';
        %         Dati.basis = 'modal';
        %
        %         X = ['    I am generating the quadrilateral mesh on L-shaped domain...']; disp(X);
        %         Tau_coarse = quad_grid_on_L_shaped_domain(2);
        %         [neighbour] = neighbours(Tau_coarse);
        %         [error] = Plot_PolyMesh(Tau_coarse,1,'k',neighbour,'-','r');
        %
        %         X = ['    I am generating the fine mesh on the same domain...']; disp(X);
        %         [Tau_fine] = quad_grid_on_L_shaped_domain(6);
        %         [neighbour] = neighbours(Tau_fine);
        %         [error] = Plot_PolyMesh(Tau_fine,2,'k',neighbour,'-','r');
        
        
        
        % L-shaped domain, fine is Voronoi, coarse is obtained by
        % agglomeration
        %         addpath(genpath('/Users/mac/Desktop/phd/POLYHEDRA/AGGLOMERATE'));
        %         addpath(genpath('/Users/mac/Desktop/phd/POLYHEDRA/AGGLOMERATE/metis-5.0.2/metismex'));
        %         melem = 1000;
        %         X = ['    I am generating the mesh with ',num2str(melem),...
        %             ' polygonal elements on L-shaped domain...']; disp(X);
        %         Tau_fine = generate_mesh_L_shaped(@MbbDomain_L_shaped,melem,100);
        %         [neighbour] = neighbours(Tau_fine);
        %         [error] = Plot_PolyMesh(Tau_fine,1,'k',neighbour,'-','r');
        %         Tau_grids{1} = Tau_fine;
        %         n_elements_coarse = Tau_fine.ne; fig_id = 1;
        %         for i=1:2
        %             X = ['    agglomeration...']; disp(X);
        %             [Tau_grids{i+1}] = agglomeration_metis(Tau_grids{i},neighbour,n_elements_coarse);
        %             [neighbour] = neighbours(Tau_grids{i+1});
        %             ne_Tau = Tau_grids{i+1}.ne;
        %             X = ['       Mesh number ',num2str(i),', num of element: ',num2str(Tau_grids{i}.ne),', num of agglomerated polygons: ',num2str(ne_Tau)];
        %             disp(X);
        %             fig_id = fig_id + 1;
        %             [error] = Plot_PolyMesh(Tau_grids{i+1},fig_id,'k',neighbour,'-','r');
        %             n_elements_coarse = Tau_grids{i+1}.ne;
        %             %Tau_fine_here = Tau_coarse;
        %         end
        %         Tau_coarse = Tau_grids{end};
        %         Tau_fine = Tau_grids{1};
        
        
        load('Tau_L_no_convex_coarse_2000_10');
        Tau_coarse = Taus{1}; Tau_fine = Taus{2};
        
        Dati.permeability = 'PolygonalChessboard';
        Dati.DiffusionCoefficientRegion = Tau_coarse; % Chessboard coarse
        %Dati.DiffusionCoefficientRegion = Tau_fine; % Chessboard fine
        Taus{1} = Tau_coarse; range_h(1) = Taus{1}.ne;% Tau_H
        Taus{2} = Tau_fine; range_h(2) = Taus{2}.ne;% Tau_h
        a=1;
        
        
        %         load('Tau_polygonal');
        %         range_h(1) = Taus{1}.ne;% Tau_H
        %         range_h(2) = Taus{2}.ne;% Tau_h
        %
        %         Dati.DiffusionCoefficientRegion = Taus{2}; % Chessboard fine
        %
        
        
    case -1
        
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
        
    case -2
        
        X = ['  The finest one is polygonal, the coarser ones are structured nested triangular.']; disp(X);
        
    case -3
        
        X = ['  The finest one is polygonal, the coarser ones are structured nested quadrilateral.']; disp(X);
        
    case -4
        
        X = ['  The finest one is polygonal, the coarser ones are triangular.']; disp(X);
        
    case -5
        
        % Try with Tau_region_airfoil_2.mat,...,Tau_region_airfoil_6.mat
        %filename_airfoil = '/Users/mac/Desktop/ALL_PHD_MATLAB_CODES/Airfoil/distmesh_WithMetis_MOX30/Tau_region_airfoil_5.mat';
        
        % Try with Airfoil_h_0p25.mat, Airfoil_h_0p125.mat or Airfoil_h_0p0625.mat
        X = ['  Airfoil profile.']; disp(X);
        filename_airfoil = ['/u/pennesi/Desktop/Multilevel_NonNested_QuadratureFree/Airfoil_poly_grids/',mesh];
        load(filename_airfoil);
        start = 2;
        %range_h = [Tau_poly{4}.ne Tau_poly{3}.ne Tau_poly{2}.ne Tau_poly{1}.ne];
        %Taus{1} = Tau_poly{4};
        %Taus{2} = Tau_poly{3};
        %Taus{3} = Tau_poly{2};
        %Taus{4} = Tau_poly{1};
        
        range_h = [Tau_poly{3}.ne Tau_poly{2}.ne];
        Taus{1} = Tau_poly{3};
        Taus{2} = Tau_poly{2};
        
    case -6
        
        % Employed only to show the p^2 rete for the non-nested case
        range_h = [512 1024];
        %range_h = [1024 2048]; start = 2;
        
        idx = 1;
        n_set = 1;%range_n_set(idx);
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
        Taus{1} = region_vec{4}; range_h(1) = Taus{1}.ne;
        
        idx = 2;
        n_set = 1;%range_n_set(idx);
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
        Taus{2} = region_vec{3}; 
        range_h(2) = Taus{2}.ne;
        
    case -7
        
        range_h = [17 1024];
        
        melem = range_h(1);
        X = ['    I am generating the mesh with ',num2str(melem),...
            ' polygonal elements...']; disp(X);
        Taus{1} = generate_mesh(@MbbDomain,melem,100);
        
        
        idx = 2;
        n_set = 1;%range_n_set(idx);
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
        Taus{2} = region_vec{2}; 
        range_h(2) = Taus{2}.ne;
        
        %save('Taus_262_18',Taus);
        %load('Taus_262_18');range_h(2) = Taus{2}.ne; range_h(1) = Taus{1}.ne;
        
    otherwise
        
        % Polygonal non-nested Voronoi grids
        
        %addpath(genpath('/u/pennesi/Desktop/PolygonClipper/'));
        for k = 1:length(range_h)
            melem = range_h(k);
            X = ['    I am generating the mesh with ',num2str(melem),...
                ' polygonal elements...']; disp(X);
            Taus{k} = generate_mesh(@MbbDomain,melem,100);
            %Taus{k} = generate_mesh(@MbbDomain_circular_crown,melem,100);
        end
        
end

for rho_black = range_rho
    for p = range_p
        
        %Dati are the same for the 2 regions
        %Dati = dati_2;
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
        %fprintf('\n');
        
        X = ['jump of rho = ',num2str(rho_black)];
        disp(X);
        fprintf('\n');
        
        for idx = start:length(range_h)
            
            %print info
            fprintf('\n');
            X = ['  Cardinality of Tau_h (finest mesh): ',num2str(range_h(idx))];
            disp(X);
            
            range = range_h(1:idx); range = range(end:-1:1);
            EndOfLoop = idx;
            
            switch nested
                
                case 1
                    
                    n_set = 1;%range_n_set(idx);
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
                    
                    for k = 1:n_file
                        range(k) = region_vec{k}.ne;
                        Taus{k} = region_vec{n_file-k+1};
                    end
                    clear region_vec; clear neigh_temp;
                    
                    EndOfLoop = 5;
                    
                case 2
                    
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
                    
                    for k = 1:n_file
                        range(k) = region_vec{k}.ne;
                        Taus{k} = region_vec{n_file-k+1};
                        %Taus{k} = region_vec{k};
                    end
                    clear region_vec; clear neigh_temp;
                    
                    EndOfLoop = 5;
                    
                case -2
                    
                    help_range = range(idx:-1:1);
                    X = ['    I am generating the mesh with ',num2str(help_range(idx)),...
                        ' elements ...']; disp(X);
                    Taus{idx} = generate_mesh(@MbbDomain,help_range(idx),100);
                    Dati.type_mesh = 'TRIA_S'; Dati.method = 'LDG'; Dati.basis = 'modal';
                    for k = length(help_range)-1:-1:1
                        melem = 2*4^(k+1);
                        X = ['    I am generating the triangular mesh with ',num2str(melem),...
                            ' elements...']; disp(X);
                        [Taus{k}] = generate_mesh_2(Dati,k+1);
                    end
                    range = 2*range; range(1) = range(1)/2;
                    EndOfLoop = idx;
                    
                case -3
                    
                    help_range = range(idx:-1:1);
                    X = ['    I am generating the mesh with ',num2str(help_range(idx)),...
                        ' elements ...']; disp(X);
                    Taus{idx} = generate_mesh(@MbbDomain,help_range(idx),100);
                    if range_h(end-1) == 4
                        index = -1;
                    elseif range_h(end-1) == 16
                        index = 0;
                    elseif range_h(end-1) == 64
                        index = 1;
                    elseif range_h(end-1) == 256
                        index = 2;
                    elseif range_h(end-1) == 1024
                        index = 3;
                    elseif range_h(end-1) == 4096
                        index = 4;
                    end
                    Dati.type_mesh = 'QUAD'; Dati.method = 'LDG'; Dati.basis = 'modal';
                    for k = length(range_h):-1:2%length(help_range)-1:-1:1
                        melem = 4^(k+index);%4^(k+1);
                        X = ['    I am generating the quadrilateral mesh with ',num2str(melem),...
                            ' elements...']; disp(X);
                        [Taus{k-1}] = generate_mesh_2(Dati,k+index);%generate_mesh_2(Dati,k+1);
                    end
                    EndOfLoop = idx;
                    
                case -4
                    
                    help_range = range(idx:-1:1);
                    X = ['    I am generating the mesh with ',num2str(help_range(idx)),...
                        ' elements ...']; disp(X);
                    Taus{idx} = generate_mesh(@MbbDomain,help_range(idx),100);
                    hmax = Taus{idx}.h*2;
                    for k = (length(help_range)-1):-1:1
                        X = ['    I am generating the triangular mesh with hmax=',num2str(hmax),...
                            '...']; disp(X);
                        [Taus{k}] = generate_UnstructTriMeshOnSquare(Dati,hmax);
                        help_range(k) = Taus{k}.ne;
                        hmax = hmax*2;
                    end
                    range = help_range(end:-1:1);
                    clear help_range;
                    EndOfLoop = idx;
                    
                otherwise
                    
                    EndOfLoop = idx;
                    
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%
            % GENERATE MATRICES %
            %%%%%%%%%%%%%%%%%%%%%
            for k = 1:EndOfLoop
                
                melem = range(k);
                
                %print info
                X = ['    Matrices for the space with ',num2str(melem),...
                    ' elements...'];
                disp(X);
                
                femregion{k} = create_dof(Dati,Taus{EndOfLoop-k+1});
                
                if k==1
                    
                    rho = get_diffusion_coefficient(femregion{k},Dati,1,rho_black);
                    femregion{k}.rho = rho;
                    
                    neighbour = neighbours(Taus{EndOfLoop-k+1});
                    %Matrices = QuadFreeMatrix2D_2(femregion{k},neighbour,P,dP,Dati);
                    Matrices = QuadFreeMatrix2D_3(femregion{k},neighbour,P,dP,Dati);
                    %Matrices = matrix2D(femregion{k},neighbour,Dati);
                    
                    rhs{k} = Matrices.f;
                    M{k} = Matrices.M;
                    A{k} = Matrices.A;
                    
                    V{k} = Matrices.V;
                    S{k} = Matrices.S;
                    I{k} = Matrices.I;
                end
                
                if (k>1) && (choice_of_A == 1)
                    neighbour = neighbours(Taus{EndOfLoop-k+1});
                    Matrices = QuadFreeMatrix2D_2(femregion{k},neighbour,P,dP,Dati);
                    %Matrices = matrix2D(femregion{k},neighbour,Dati);
                    
                    rhs{k} = Matrices.f;
                    M{k} = Matrices.M;
                    A{k} = Matrices.A;
                    
                    V{k} = Matrices.V;
                    S{k} = Matrices.S;
                    I{k} = Matrices.I;
                end
                
            end
            clear Matrices, clear neighbour;
            f = rhs{1}; clear rhs;
            
            %Generate projectors and matrix for sublevels and solve with Additive Schwarz PCG
            M_inv = inv_SIPDG_MassMatrix(M{1}, ...
                femregion{1}.ne, femregion{1}.nln);
            switch choice_of_size
                
                case 1
                    
                    % h < H < HH
                    for k = 2:EndOfLoop-1
                        
                        local_solvers = MY_LOCAL_SOLVERS_2(femregion{k}, femregion{1}, A{1});
                        M_hH = Mixed_MassMatrix_QuadFree(femregion{k+1}, femregion{1}, P, p);
                        
                        %Prolongation operator R_h{k-1} : femregion{k} --> femregion{k-1}
                        R_0 = M_inv * M_hH;
                        
                        if (choice_of_A == 2)
                            V{k+1} = (R_0')*V{1}*R_0;
                            S{k+1} = (Taus{EndOfLoop}.h/Taus{EndOfLoop-k}.h)*(R_0')*S{1}*R_0;
                            I{k+1} = (R_0')*I{1}*R_0;
                            A{k+1} = V{k+1} - transpose(I{k+1}) - I{k+1} + S{k+1};
                        elseif (choice_of_A == 3)
                            A{k+1} = (R_0')*A{1}*R_0;
                        end
                        
                        coarse_solver = MY_COARSE_SOLVER(femregion{k+1}, A{k+1}, R_0, femregion{k}.ne);
                        if p==range_p(1)
                            N_h(EndOfLoop,k) = range(1);
                            N_H(EndOfLoop,k) = range(k);
                            N_HH(EndOfLoop,k) = range(k+1);
                        end
                        
                        % Solve with Additive Schwarz PCG
                        z0 = zeros(femregion{1}.ndof,1);
                        %tstart = tic;
                        X = ['    solving...']; disp(X);
                        [uh,~,iter_pcg(p,EndOfLoop,k),~,condi(p,EndOfLoop,k),~] = my_pcg(A{1}, z0, f, 100000, tol, local_solvers, coarse_solver);
                        %telapse = toc(tstart);
                        %iter = iter_pcg(p,EndOfLoop,k);
                        %X = ['     - ASPCG completed in ',...
                        %    num2str(telapse),' for N_h = ',num2str(range(1)),', N_H = ',num2str(range(k)),', with iter=',num2str(iter)];
                        %disp(X);
                        X = ['     - N_h = ',num2str(N_h(EndOfLoop,k)),...
                            ', N_H = ',num2str(N_H(EndOfLoop,k)),...
                            ', N_HH = ',num2str(N_HH(EndOfLoop,k)),...
                            ', K(P_ad) = ',num2str(condi(p,EndOfLoop,k)),...
                            ' it = ',num2str(iter_pcg(p,EndOfLoop,k))];
                        disp(X);
                        
                    end
                    
                case 2
                    
                    % h = H < HH
                    local_solvers = MY_LOCAL_SOLVERS(femregion{1}, A{1});
                    for k = startLine615:EndOfLoop
                        
                        %M_hH = Mixed_MassMatrix_QuadFree(femregion{k}, femregion{1}, P, p);
                        
                        if q == 0
                            q = p;
                        end
                        nln=0.5.*(q+1).*(q+2);
                        femregion{k}.fem = q;
                        femregion{k}.nln = nln;
                        femregion{k}.ndof = nln*femregion{k}.ne;
                        M_hH = Mixed_MassMatrix_QuadFree_2(femregion{k}, femregion{1});
                        
                        
                        %Prolongation operator R_h{k-1} : femregion{k} --> femregion{k-1}
                        R_0 = M_inv * M_hH;
                        
                        if (choice_of_A == 2)
                            V{k} = (R_0')*V{1}*R_0;
                            S{k} = (Taus{EndOfLoop}.h/Taus{EndOfLoop-k+1}.h)*(R_0')*S{1}*R_0;
                            I{k} = (R_0')*I{1}*R_0;
                            A{k} = V{k} - transpose(I{k}) - I{k} + S{k};
                        elseif (choice_of_A == 3)
                            A{k} = (R_0')*A{1}*R_0;
                        end
                        
                        coarse_solver = MY_COARSE_SOLVER(femregion{k}, A{k}, R_0, femregion{1}.ne);
                        if p==range_p(1)
                            N_h(EndOfLoop,k) = range(1);
                            N_H(EndOfLoop,k) = range(1);
                            N_HH(EndOfLoop,k) = range(k);
                        end
                        
                        % Solve with Additive Schwarz PCG
                        z0 = zeros(femregion{1}.ndof,1);
                        %tstart = tic;
                        X = ['    solving...']; disp(X);
                        [uh,~,iter_pcg(p,EndOfLoop,k),~,condi(p,EndOfLoop,k),~] = my_pcg(A{1}, z0, f, 100000, tol, local_solvers, coarse_solver);
                        %telapse = toc(tstart);
                        %iter = iter_pcg(p,EndOfLoop,k);
                        %X = ['     - ASPCG completed in ',...
                        %    num2str(telapse),' for N_h = ',num2str(range(1)),', N_H = ',num2str(range(k)),', with iter=',num2str(iter)];
                        %disp(X);
                        X = ['     - N_h = ',num2str(N_h(EndOfLoop,k)),...
                            ', N_H = ',num2str(N_H(EndOfLoop,k)),...
                            ', N_HH = ',num2str(N_HH(EndOfLoop,k)),...
                            ', K(P_ad) = ',num2str(condi(p,EndOfLoop,k)),...
                            ' it = ',num2str(iter_pcg(p,EndOfLoop,k))];
                        disp(X);
                        
                    end
                    
                case 3
                    
                    % h < H = HH
                    for k = 2:EndOfLoop
                        
                        local_solvers = MY_LOCAL_SOLVERS_2(femregion{k}, femregion{1}, A{1});
                        M_hH = Mixed_MassMatrix_QuadFree(femregion{k}, femregion{1}, P, p);
                        
                        %Prolongation operator R_h{k-1} : femregion{k} --> femregion{k-1}
                        R_0 = M_inv * M_hH;
                        
                        if (choice_of_A == 2)
                            V{k} = (R_0')*V{1}*R_0;
                            S{k} = (Taus{EndOfLoop}.h/Taus{EndOfLoop-k+1}.h)*(R_0')*S{1}*R_0;
                            I{k} = (R_0')*I{1}*R_0;
                            A{k} = V{k} - transpose(I{k}) - I{k} + S{k};
                        elseif (choice_of_A == 3)
                            A{k} = (R_0')*A{1}*R_0;
                        end
                        
                        coarse_solver = MY_COARSE_SOLVER(femregion{k}, A{k}, R_0, femregion{k}.ne);
                        if p==range_p(1)
                            N_h(EndOfLoop,k) = range(1);
                            N_H(EndOfLoop,k) = range(k);
                            N_HH(EndOfLoop,k) = range(k);
                        end
                        
                        % Solve with Additive Schwarz PCG
                        z0 = zeros(femregion{1}.ndof,1);
                        %tstart = tic;
                        X = ['    solving...']; disp(X);
                        [uh,~,iter_pcg(p,EndOfLoop,k),~,condi(p,EndOfLoop,k),~] = my_pcg(A{1}, z0, f, 100000, tol, local_solvers, coarse_solver);
                        %telapse = toc(tstart);
                        %iter = iter_pcg(p,EndOfLoop,k);
                        %X = ['     - ASPCG completed in ',...
                        %    num2str(telapse),' for N_h = ',num2str(range(1)),', N_H = ',num2str(range(k)),', with iter=',num2str(iter)];
                        %disp(X);
                        X = ['     - N_h = ',num2str(N_h(EndOfLoop,k)),...
                            ', N_H = ',num2str(N_H(EndOfLoop,k)),...
                            ', N_HH = ',num2str(N_HH(EndOfLoop,k)),...
                            ', K(P_ad) = ',num2str(condi(p,EndOfLoop,k)),...
                            ' it = ',num2str(iter_pcg(p,EndOfLoop,k))];
                        disp(X);
                        
                    end
                    
            end
            clear femregion; clear Tau; clear A; ...
                clear V; clear I; clear S; clear M;
        end
    end
end % loop over rho_black
% Save info:
path = 'AS_RESULTS/';
if(nested >= 0)
    filename1 = ['Nest',num2str(nested)];
elseif(nested == -5)
    filename1 = ['Airfoil_',filename_airfoil(83:end)];
else
    filename1 = ['Nest_m_',num2str(abs(nested))];
end

%if nested == -6
%    filename1 = [filename1,'_2'];
%end

filename1 = [filename1,'_',num2str(Taus{2}.ne),'_q',num2str(q)];

file = [char(path),char(filename1)];
save(file,'iter_pcg','condi','N_h','N_H','N_HH');
