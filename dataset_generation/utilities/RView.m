classdef RView < handle & Utility
    properties
        mesh
        loc2glob
        glob2loc
        elem_num
    end
    methods
        function view = RView(mesh,ids,dim)
            if nargin == 1
                dim = manifold(mesh);
                ids = 1:mesh.elem_num(dim+1);
                if isempty(dim) % => isempty(ids)
                    dim = 0;
                end
            end

            view.mesh = mesh;

            % loc2glob
            view.loc2glob = cell(1,mesh.dim+1);
            view.loc2glob{dim+1} = ids;
            if dim > 0
                for d = dim:-1:2
                    view.loc2glob{d} = unique(get2(mesh,view.loc2glob{d+1},d));
                end
                view.loc2glob{1} = reshape(unique(get1(mesh,view.loc2glob{2})),1,[]);
            end

            % glob2loc
            view.glob2loc = cell(1,mesh.dim+1);
            for i = 1:dim+1
                view.glob2loc{i} = sparse([]);
                view.glob2loc{i}(view.loc2glob{i}) = 1:length(view.loc2glob{i});
            end

            % elem_num
            view.elem_num = zeros(1,mesh.dim+1);
            for i = 1:dim+1
                view.elem_num(i) = length(view.loc2glob{i});
            end
        end
        function t = tol(view)
           t = view.mesh.tol;
        end % can get but not set
        function d = dim(view)
           d = view.mesh.dim;
        end % can get but not set
        function c = convex(view)
           c = view.mesh.convex;
        end % can get but not set
        function id = add2view(view,mesh_id,dim)
            % does NOT work recursively:
            % use merge in that case
            view.elem_num(dim+1) = view.elem_num(dim+1) + 1;
            id = view.elem_num(dim+1);
            view.loc2glob{dim+1}(id) = mesh_id;
            view.glob2loc{dim+1}(mesh_id) = id;
        end
        function ids = get_share(view,id,dim)
            mesh_ids = get_share(view.mesh,view.loc2glob{dim+1}(id),dim);
            % return elements inside view only
            ids = view.glob2loc{dim+2}(mesh_ids(ismember(mesh_ids,view.loc2glob{dim+2})));
        end
        function mesh = copy(view)
            mesh = RMesh(view.dim);
            mesh.convex = view.convex;
            for i = 1:view.elem_num(1)
                add_new0(mesh,get0(view,i));
            end
            for i = 1:view.elem_num(2)
                add_new1(mesh,get1(view,i));
            end
            for d = 2:view.dim
                for i = 1:view.elem_num(d+1)
                    add_new2(mesh,get2(view,i,d),d);
                end
            end
        end
        function id = add_poly2D(view,poly)
            n = length(poly);
            elem = zeros(1,n);
            for i = 1:n
                i_plus = mod(i,n)+1;
                edge = [poly(i),poly(i_plus)];
                id = add1(view,edge);
                elem(i) = id;
            end
            id = add2(view,elem,2);
        end
        function to_view = merge(to_view,from_mesh)
            if isempty(manifold(from_mesh))
                return
            end

            % change vert dimension
            vert = vertices(from_mesh);
            md = min(from_mesh.dim,to_view.dim);
            Md = max(from_mesh.dim,to_view.dim);
            vert = [vert(:,1:md),zeros(from_mesh.elem_num(1),Md-md)];

            new_ids = zeros(1,from_mesh.elem_num(1));
            for i = 1:from_mesh.elem_num(1)
                new_ids(i) = add0(to_view,vert(i,:));
            end

            old_ids = new_ids;

            egs = edges(from_mesh);
            egs = old_ids(egs);
            new_ids = zeros(1,from_mesh.elem_num(2));
            for i = 1:from_mesh.elem_num(2)
                new_ids(i) = add1(to_view,egs(i,:));
            end


            for d = 2:from_mesh.dim
                old_ids = new_ids;
                new_ids = zeros(1,from_mesh.elem_num(d+1));
                for i = 1:from_mesh.elem_num(d+1)
                    elem = get2(from_mesh,i,d);
                    elem = old_ids(elem);
                    new_ids(i) = add2(to_view,elem,d);
                end
            end

        end
        function pyramid(base_view,point)
        % base needs to have the same dimension of pyramid!
        dim = manifold(base_view);
        if dim == 0
            % in this case base_mesh is made only of 1 point
            id = add0(base_view,point);
            add1(base_view,[1 id]);
            return
        end

        collinear = @(vert,point) rank(vert-point,base_view.dim*base_view.tol) <= dim;
        bndry = boundary(base_view);

        if collinear(vertices(bndry),point)
            elem = 1:bndry.elem_num(dim);
            if dim == 1
                add1(bndry,elem);
            else
                add2(bndry,elem,dim);
            end
            add2(base_view,1:base_view.elem_num(dim+1),dim+1);
            return
        end

        graph = connectivity(bndry,dim-1);
        coll_edges = false(1,graph.elem_num(2));

        for i = 1:graph.elem_num(2)
            edge = get1(graph,i);
            view = RView(bndry,edge,dim-1);
            if collinear(vertices(view),point)
               coll_edges(i) = true;
            end
        end

        clust_graph = RView(graph,find(coll_edges),1);
        [clust,n_clust] = clusters(clust_graph);

        for i = 1:n_clust
            I = clust_graph.loc2glob{1}(clust == i);
            view = RView(bndry,I,dim-1);
            pyramid(view,point);
        end

        for i = setdiff(1:graph.elem_num(1),clust_graph.loc2glob{1})
            view = RView(bndry,i,dim-1);
            pyramid(view,point);
        end

        add2(base_view,1:base_view.elem_num(dim+1),dim+1);

        end


        % set methods
        function set0(view,points,ids)
            set0(view.mesh,points,view.loc2glob{1}(ids));
        end
        function set1(view,edges,ids)
            set1(view.mesh,edges,view.loc2glob{2}(ids));
        end
        function set2(view,elem,id,dim)
            set2(view.mesh,elem,view.loc2glob{dim+1}(id),dim);
        end

        % get methods
        function elem = get0(view,ids)
            % since I am getting the points coordinates
            % there is no global to local indexing
            elem = get0(view.mesh,view.loc2glob{1}(ids));
        end
        function elem = get1(view,ids)
            mesh_elem = get1(view.mesh,view.loc2glob{2}(ids));
            elem = view.glob2loc{1}(mesh_elem);
            % elem = full(view.glob2loc{1}(mesh_elem));
        end
        function elem = get2(view,ids,dim)
            mesh_elem = get2(view.mesh,view.loc2glob{dim+1}(ids),dim);
            elem = view.glob2loc{dim}(mesh_elem);
            % elem = full(view.glob2loc{dim}(mesh_elem));
        end

        % add methods
        function id = add_new0(view,point)
            mesh_id = add_new0(view.mesh,point);
            id = add2view(view,mesh_id,0);
        end
        function id = add_new1(view,edge)
            mesh_id = add_new1(view.mesh,view.loc2glob{1}(edge));
            id = add2view(view,mesh_id,1);
        end
        function id = add_new2(view,elem,dim)
            mesh_id = add_new2(view.mesh,view.loc2glob{dim}(elem),dim);
            id = add2view(view,mesh_id,dim);
        end
        function [id,TF] = add0(view,point)
            id = 1;
            while id <= view.elem_num(1) && ...
                    norm(point - get0(view,id),inf) > view.tol
                id = id + 1;
            end
            TF = id > view.elem_num(1);
            if TF
                mesh_id = add0(view.mesh,point);
                add2view(view,mesh_id,0);
            end
        end
        function [id,TF] = add1(view,edge)
            edge = sort(edge);
            id = 1;
            while id <= view.elem_num(2) && ...
                    not(isequal(edge,get1(view,id)))
                id = id + 1;
            end
            TF = id > view.elem_num(2);
            if TF
                mesh_id = add1(view.mesh,view.loc2glob{1}(edge));
                add2view(view,mesh_id,1);
            end
        end
        function [id,TF] = add2(view,elem,dim)
            elem = sort(elem);
            id = 1;
            while id <= view.elem_num(dim+1) && ...
                    not(isequal(elem,get2(view,id,dim)))
                id = id + 1;
            end
            TF = id > view.elem_num(dim+1);
            if TF
                mesh_id = add2(view.mesh,view.loc2glob{dim}(elem),dim);
                add2view(view,mesh_id,dim);
            end
        end

        % replace methods
        function id_remove = replace1(view,id,add_fun,N_add)
            mesh_id_remove = replace1 ...
                (view.mesh,view.loc2glob{2}(id),add_fun,N_add);

            id_remove = view.glob2loc{2}(mesh_id_remove);
            view.loc2glob{2}(id_remove) = [];
            view.glob2loc{2}(mesh_id_remove) = 0; % glob2loc is sparse
            view.elem_num(2) = view.elem_num(2) - 1;
        end
        function id_remove = replace2(view,id,add_fun,N_add,dim)
            mesh_id_remove = replace2 ...
                (view.mesh,view.loc2glob{dim+1}(id),add_fun,N_add,dim);

            id_remove = view.glob2loc{dim+1}(mesh_id_remove);
            view.loc2glob{dim+1}(id_remove) = [];
            view.glob2loc{dim+1}(mesh_id_remove) = 0; % glob2loc is sparse
            view.elem_num(dim+1) = view.elem_num(dim+1) - 1;
        end

    end

end
