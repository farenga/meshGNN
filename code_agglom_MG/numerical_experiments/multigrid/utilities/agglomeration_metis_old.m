function [region_coarse] = agglomeration_metis_old(region,neighbour,n_elements_coarse)
    
    edges = region.edges_phys;
    nodes = region.coord;
    
    % METIS CALL
    eind = [];
    eptr(1) = 0;
    A = sparse(region.ne,region.ne);
    for i = 1:region.ne
        elements = neighbour.neigh{i};
        elements(elements<0) = [];
        A(i,elements) = 1;
%         eptr(i+1) = eptr(i)+length(region.connectivity{i});
%         eind = [eind region.connectivity{i}-1];
%         vsize(i) = polyarea(region.coords_element{i}(:,1),region.coords_element{i}(:,2));
    end
    epart = metismex('PartGraphKway',A,n_elements_coarse);epart = epart+1;
    t = sort(unique(epart));
    t_check = [1:max(epart)];
    difference = setdiff(t_check,t);
    
    if isempty(difference)==false
        for i = 1:length(difference)
            epart(epart>difference(i)-i+1) = epart(epart>difference(i)-i+1)-1;
        end
    end
    
    nparts = max(epart);
    region_temp.elements= cell(nparts,1);  % fine elements belonging to each part
    for i = 1:length(epart)
        region_temp.elements{epart(i)} = union(region_temp.elements{epart(i)},i);
    end   
    
    region_temp.connectivity= cell(nparts,1);
    region_temp.edges= cell(nparts,1);
    for i = 1:nparts
        for j = 1:length(region_temp.elements{i})
            check_neigh = 0;
            neigh = neighbour.neigh{region_temp.elements{i}(j)};
            for k = 1:length(neigh)
                if neigh(k)==-1 || i~=epart(neigh(k))
                    region_temp.connectivity{i} = union(region_temp.connectivity{i},edges{region_temp.elements{i}(j)}(k,:));
                    region_temp.edges{i}=[ region_temp.edges{i};edges{region_temp.elements{i}(j)}(k,:)];
                else
                    check_neigh = 1;
                end
            end
            if check_neigh == 0
                elem = region_temp.elements{i};
                figure
                for s = elem
                    coordinates = region.coords_element{s};
                    patch(coordinates(:,1),coordinates(:,2),'w');
                    hold on
                    axis equal
                end
                error('Disconnected elements given by METIS, Agglomerated element: %d',i);
                
            end
        end
        region_temp.connectivity{i} = unique(region_temp.connectivity{i});
        
    end        

    % Counterclockwise nodes order
    for i = 1:nparts
        nodes_new = region_temp.edges{i}(1,:);
        set = nodes_new(end);
        while (length(nodes_new)<length(region_temp.connectivity{i}))
            index = find(region_temp.edges{i}(:,1)==set);
            nodes_new = [nodes_new region_temp.edges{i}(index,2)];
            set = nodes_new(end);
        end
        region_temp.connectivity{i} = nodes_new;
    end
    
    % Check collinearity
    nodes_ID = [];    
    for i = 1:size(region_temp.connectivity,1)
        ne_i = length(region_temp.connectivity{i});
        flag = ones(ne_i,1);
        for j = 1:ne_i
            if j<ne_i-1
                v1 = nodes(region_temp.connectivity{i}(j),:);
                v2 = nodes(region_temp.connectivity{i}(j+1),:);
                v3 = nodes(region_temp.connectivity{i}(j+2),:);
                local_ID = [j j+1 j+2];
            elseif j==ne_i-1
                v1 = nodes(region_temp.connectivity{i}(j),:);
                v2 = nodes(region_temp.connectivity{i}(j+1),:);
                v3 = nodes(region_temp.connectivity{i}(1),:);
                local_ID = [j j+1 1];
            elseif j==ne_i
                v1 = nodes(region_temp.connectivity{i}(j),:);
                v2 = nodes(region_temp.connectivity{i}(1),:);
                v3 = nodes(region_temp.connectivity{i}(2),:);
                local_ID = [j 1 2];
            end
            if collinear(v1,v2,v3)==1
                flag(local_ID(2))=0;
            end
                
        end
        region_temp.connectivity{i}(flag==0)=[];
        nodes_ID = union(nodes_ID, region_temp.connectivity{i});
    end
        
    nodes_ID = unique(nodes_ID);
        
    % Renumbering
    ID_map = zeros(1,size(nodes,1));
    for i = 1:length(nodes_ID)
        ID_map(nodes_ID(i)) = i; 
    end    
    region_coarse.connectivity= cell(nparts,1);
    region_coarse.edges= cell(nparts,1);
    
    for i = 1:nparts
        for j = 1:length(region_temp.connectivity{i})
            region_coarse.connectivity{i}(j) = ID_map(region_temp.connectivity{i}(j));
        end
    end
    region_coarse.coord = nodes(nodes_ID,:);
    coords = region_coarse.coord;
    
    for i = 1:length(region_coarse.connectivity)
        nedge(i) = length(region_coarse.connectivity{i});
        local_coords = coords(region_coarse.connectivity{i},:);
        x_min = min(local_coords(:,1));x_max = max(local_coords(:,1));
        y_min = min(local_coords(:,2));y_max = max(local_coords(:,2));
        BBox(i,:)=[x_min x_max y_min y_max];
        coords_element{i} = local_coords;
        [centroids(i,:)]=get_area_centroid_weights( coords_element{i});
         edges{i} = [[1:nedge(i)]', [2:nedge(i) 1]'];
         edges_phys{i} = [region_coarse.connectivity{i}' region_coarse.connectivity{i}([2:end 1])'];
    end
    Element = region_coarse.connectivity;
    
    edges_ID = cell2mat(edges_phys');
    repeated_edges = zeros(size(edges_ID,1),1);
    for i = 1:size(edges_ID,1)
        if repeated_edges(i)== 0
            v1 = edges_ID(i,1); v2 = edges_ID(i,2);
            for j = i+1:size(edges_ID,1)
                if v1 == edges_ID(j,2) && v2 == edges_ID(j,1);
                    repeated_edges(j) = 1;
                end
            end
        end
    end
    
    ne = length(Element);
    edges_ID(repeated_edges==1,:)=[];
    elem_area = zeros(ne,1);
    max_kb = cell(ne,1);
    for i = 1:ne
        elem_area(i) = polyarea(coords_element{i}(:,1),coords_element{i}(:,2));
        max_kb{i} = zeros(size(edges_phys{i},1));
        for j = 1:size(edges_phys{i},1)
            v1 = coords(edges_phys{i}(j,1),:); v2 = coords(edges_phys{i}(j,2),:);
            for k = 1:size(edges_phys{i},1)
                if k~=j
                    v3 = coords(edges_phys{i}(k,2),:);
                    [x_tria,y_tria]=poly2cw([v1(1) v2(1) v3(1)],[v1(2) v2(2) v3(2)]);
                    [x1,y1] = polybool('intersection',coords_element{i}(end:-1:1,1),coords_element{i}(end:-1:1,2),x_tria,y_tria);
                    area = polyarea(x_tria,y_tria);
                    if any(isnan(x1))== false && abs(polyarea(x1,y1)- area)<1e-13
                        if area>max_kb{i}(j)
                            max_kb{i}(j) = area;
                        end
                    end
                end
            end

        end
    end


    region_coarse.nedges= nedge';
    region_coarse.ne=length(region_coarse.connectivity);
    region_coarse.BBox = BBox;
    region_coarse.coords_element=coords_element;
    region_coarse.edges = edges;
    region_coarse.edges_phys = edges_phys;
    region_coarse.centroids = centroids;
    region_coarse.edges_ID=edges_ID;
    region_coarse.aggl_elements = region_temp.elements;
    region_coarse.area = elem_area;
    region_coarse.max_kb = max_kb;