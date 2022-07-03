clear all, close all

id_fig = 0;
for n_set = 1:5
    
    Dati = dati;
    n_file = 5;
    set_str = strcat('_MGG_set',num2str(n_set));
    
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
    
    for k = 1:5
        id_fig = id_fig + 1;
        region_vec{k}.ne
        [error] = Plot_PolyMesh(region_vec{k},id_fig,'k',neigh_vec{k},'b');
    end
end

%for k = 1:5
%    region_vec{k}.ne
%    [error] = Plot_PolyMesh(region_vec{k},k,'k',neigh_vec{k},'b');
%end
%%

%==========================================
%  Main of AS preconditioner for CG method
%==========================================
clear all;
close all;

id_fig = 0;
tol = 1e-8;
%smooth = [3 5 8 12 16 20];

range_h = [16 64 256 1024]; start = 4;

nested = -4;

% Generate all the grids needed for the numerical tests
if(nested==1)
    
    %range_h = [512 1024 2048 4096]; start = 1;
    %range_n_set = [1 1 1 1];
    range_h = [2048]; start = 1;
    range_n_set = [1];
    
    
elseif nested == 2
    
    Dati = dati;
    Dati.type_mesh = 'QUAD';
    Dati.method = 'LDG';
    Dati.basis = 'modal';
    
    for k = 2:length(range_h)+1
        melem = 4^k;
        X = ['    I am building the quadrilateral mesh with ',num2str(melem),...
            ' elements...']; disp(X);
        [Taus{k-1}] = generate_mesh_2(Dati,k);
        [neighbor] = neighbours(Taus{k-1});
        id_fig = id_fig + 1;
        [error] = Plot_PolyMesh(Taus{k-1},id_fig,'k',neighbor,'k');
    end
    
elseif nested == 3
    
    Dati = dati;
    Dati.type_mesh = 'TRIA_S';
    Dati.method = 'LDG';
    Dati.basis = 'modal';
    range_h = 2*range_h;
    
    for k = 2:length(range_h)+1
        melem = 2^(2*k+1);
        X = ['    I am generating the triangular mesh with ',num2str(melem),...
            ' elements...']; disp(X);
        [Taus{k-1}] = generate_mesh_2(Dati,k);
        [neighbor] = neighbours(Taus{k-1});
        id_fig = id_fig + 1;
        [error] = Plot_PolyMesh(Taus{k-1},id_fig,'k',neighbor,'k');
    end
    
elseif nested == 4
    
    % unstructured triangles
    %range_h = [582 1086 2198 4318]; start = 1;
    range_h = [582]; start = 1;
    
    % for each fine mesh you have 4 nested grids to be loaded
    
elseif nested == -3
    
    X = ['Fine poly, coarse struct quad.']; disp(X);
    
elseif nested == -4
    
    X = ['Fine poly, coarse struct tria.']; disp(X);
    
else
    
    %addpath(genpath('/u/pennesi/Desktop/PolygonClipper/'));
    for k = 1:length(range_h)
        
        melem = range_h(k);
        
        %print info
        X = ['    I am generating the mesh with ',num2str(melem),...
            ' elements...']; disp(X);
        
        Taus{k} = generate_mesh(@MbbDomain,melem,100);
    end
end


for idx = start:length(range_h)
    
    %print info
    fprintf('\n');
    X = ['  Cardinality of Tau_h (finest mesh): ',num2str(range_h(idx))];
    disp(X);
    
    range = range_h(1:idx); range = range(end:-1:1);
    
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
        
        
        for n_reg = 1:n_file
            ne = region_vec{n_reg}.ne;
            for i = 1:ne
                coord = region_vec{n_reg}.coords_element{i};
                cc = 1;
                for j = 1:size(coord,1)-1
                    for k = j+1:size(coord,1)
                        distance(cc) = sqrt((coord(j,1)-coord(k,1)).^2+(coord(j,2)-coord(k,2)).^2);
                        cc = cc+1;
                    end
                end
                H(i) = max(distance);
            end
            
            h_max(n_reg) = max(H);
            h_min(n_reg) = min(H);
            clear H;
        end
        h_max(2:end)./h_max(1)
        
        for k = 1:n_file
            id_fig = id_fig + 1;
            region_vec{k}.ne
            %[error] = Plot_PolyMesh(region_vec{k},id_fig,'0',neigh_vec{k},'k');
            %axis off;
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
        
        % Check the value of H/h
        for n_reg = 1:n_file
            ne = region_vec{n_reg}.ne;
            for i = 1:ne
                coord = region_vec{n_reg}.coords_element{i};
                cc = 1;
                for j = 1:size(coord,1)-1
                    for k = j+1:size(coord,1)
                        distance(cc) = sqrt((coord(j,1)-coord(k,1)).^2+(coord(j,2)-coord(k,2)).^2);
                        cc = cc+1;
                    end
                end
                H(i) = max(distance);
            end
            
            h_max(n_reg) = max(H);
            h_min(n_reg) = min(H);
            clear H;
        end
        h_max(2:end)./h_max(1)
        
        for k = 1:n_file
            id_fig = id_fig + 1;
            region_vec{k}.ne
            %[error] = Plot_PolyMesh(region_vec{k},id_fig,'k',neigh_vec{k},'k');
            axis off;
        end
        clear region_vec; clear neigh_temp;
        
    elseif nested == -3
        
        help_range = range(idx:-1:1);
        X = ['    I am generating the mesh with ',num2str(help_range(idx)),...
            ' elements ...']; disp(X);
        Taus{idx} = generate_mesh(@MbbDomain,help_range(idx),100);
        Dati = dati;
        Dati.type_mesh = 'QUAD'; Dati.method = 'LDG'; Dati.basis = 'modal';
        for k = length(help_range)-1:-1:1
            melem = 4^(k+1);
            X = ['    I am generating the quadrilateral mesh with ',num2str(melem),...
                ' elements...']; disp(X);
            [Taus{k}] = generate_mesh_2(Dati,k+1);
        end
        EndOfLoop = idx;
        
        for k = 1:length(help_range)
            id_fig = id_fig + 1;
            neighbour = neighbours(Taus{k});
            [error] = Plot_PolyMesh(Taus{k},id_fig,'k',neighbour,'k');
            axis off;
        end
        
    elseif nested == -4
        
        help_range = range(idx:-1:1);
        X = ['    I am generating the mesh with ',num2str(help_range(idx)),...
            ' elements ...']; disp(X);
        Taus{idx} = generate_mesh(@MbbDomain,help_range(idx),100);
        Dati = dati;
        Dati.type_mesh = 'TRIA_S'; Dati.method = 'LDG'; Dati.basis = 'modal';
        for k = length(help_range)-1:-1:1
            melem = 2*4^(k+1);
            X = ['    I am generating the triangular mesh with ',num2str(melem),...
                ' elements...']; disp(X);
            [Taus{k}] = generate_mesh_2(Dati,k+1);
        end
        range = 2*range; range(1) = range(1)/2;
        EndOfLoop = idx;
        
        for k = 1:length(help_range)
            id_fig = id_fig + 1;
            neighbour = neighbours(Taus{k});
            [error] = Plot_PolyMesh(Taus{k},id_fig,'k',neighbour,'k');
            axis off;
        end
        
    end
    
    
end



