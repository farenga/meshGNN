classdef RMesh < matlab.mixin.Copyable & Utility
    % RMesh Refinable Mesh Class
    properties
        dim        
        convex = true
        tol = 1e-12 % infinity norm vertex tolerance
        elem_num
        elem
        share
    end
    methods
        function mesh = RMesh(dim)
            mesh.dim = dim;
            mesh.elem_num = zeros(1,dim+1);
            mesh.elem = cell(1,dim+1);
            mesh.share = cell(1,dim+1);
        end
        function ids = get_share(mesh,id,dim)
            ids = mesh.share{dim+1}{id};
        end
        
        % different versions of the same methods are used
        % in order to avoid if statements and improve efficiency
        
        % set methods
        function set0(mesh,points,ids)
            mesh.elem{1}(ids,:) = points;
        end
        function set1(mesh,edges,ids)
            mesh.elem{2}(ids,:) = sort(edges);
        end
        function set2(mesh,elem,id,dim)
            mesh.elem{dim+1}{id} = sort(elem);
        end
        
        % get methods
        function elem = get0(mesh,ids)
            elem = mesh.elem{1}(ids,:);
        end
        function elem = get1(mesh,ids)
            elem = mesh.elem{2}(ids,:);
        end
        function elem = get2(mesh,ids,dim)
            elem = [mesh.elem{dim+1}{ids}];
        end
        
        % add methods
        function id = add_new0(mesh,point)
            mesh.elem_num(1) = mesh.elem_num(1)+1;
            id = mesh.elem_num(1);
            set0(mesh,point,id);
            mesh.share{1}{id} = [];
        end
        function id = add_new1(mesh,edge)
            mesh.elem_num(2) = mesh.elem_num(2)+1;
            id = mesh.elem_num(2);
            set1(mesh,edge,id);
            mesh.share{2}{id} = [];
            for v = edge
                mesh.share{1}{v}(end+1) = id;
            end
        end
        function id = add_new2(mesh,elem,dim)
            mesh.elem_num(dim+1) = mesh.elem_num(dim+1)+1;
            id = mesh.elem_num(dim+1);
            set2(mesh,elem,id,dim);
            mesh.share{dim+1}{id} = [];
            for e = elem
                mesh.share{dim}{e}(end+1) = id;
            end
        end
        
        % replace methods
        function id_remove = replace1(mesh,id,add_fun,N_add)
        % last edges are used to replace edge get1(mesh,id)
        
            % the last element replaces the old one
            % so add new elements in this order:
            for i = 2:N_add
            	add_fun(i);
            end
            add_fun(1);
        
            % replace old edge with last edge
            old = get1(mesh,id);
            id_last = mesh.elem_num(2);
            last = get1(mesh,id_last);
            set1(mesh,last,id);
            for e = old
                mesh.share{1}{e}(mesh.share{1}{e} == id) = [];
            end
            for e = last
                mesh.share{1}{e}(mesh.share{1}{e}==id_last) = id;
            end

            % remove last edge
            id_remove = id_last;
            mesh.elem{2}(id_remove,:)= [];
            mesh.share{2}(id_last) = [];
            mesh.elem_num(2) = mesh.elem_num(2) - 1;            
            
            % replace elements used by the same using old
            new_ids = (id_last-N_add+1):(id_last-1);
            mesh.share{2}(new_ids) = mesh.share{2}(id);
            
            % elements using old elem now use replace elements
            for s = mesh.share{2}{id}
                el = get2(mesh,s,2);
                set2(mesh,[el,new_ids],s,2);
            end 
        end
        function id_remove = replace2(mesh,id,add_fun,N_add,dim)
        % last mesh elements are used to replace get2(mesh,id,dim)
            
            % the last element replaces the old one
            % so add new elements in this order:
            for i = 2:N_add
            	add_fun(i);
            end
            add_fun(1);
            
            % replace old element with the last element
            old = get2(mesh,id,dim);
            id_last = mesh.elem_num(dim+1);
            last = get2(mesh,id_last,dim);
            set2(mesh,last,id,dim);
            for e = old
                mesh.share{dim}{e}(mesh.share{dim}{e} == id) = [];
            end
            for e = last
                mesh.share{dim}{e}(mesh.share{dim}{e}==id_last) = id;
            end

            % remove last element
            id_remove = id_last;
            mesh.elem{dim+1}(id_remove)= [];
            mesh.share{dim+1}(id_last) = [];
            mesh.elem_num(dim+1) = mesh.elem_num(dim+1) - 1;

            % replace elements used by the same using old
            new_ids = (id_last-N_add+1):(id_last-1);
            mesh.share{dim+1}(new_ids) = mesh.share{dim+1}(id);
            
            % elements using old elem now use replace elements
            for s = mesh.share{dim+1}{id}
                el = get2(mesh,s,dim+1);
                set2(mesh,[el,new_ids],s,dim+1);
            end 
        end

    end
    methods (Access = ?RView)         
        function id = add0(mesh,point)
        	id = add_new0(mesh,point); 
        end
        function id = add1(mesh,edge)
            id = add_new1(mesh,edge);
        end
        function id = add2(mesh,elem,dim)
        	id = add_new2(mesh,elem,dim);
        end
        
    end
end