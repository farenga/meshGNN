classdef (Abstract) Utility < handle
    methods
        function move(mesh,v)
            new_vert = vertices(mesh)+v;
            vertices(mesh,new_vert);
        end
        function scale(mesh,h)
            new_vert = vertices(mesh)*h;
            vertices(mesh,new_vert);
        end
        function normalize(mesh)
            v = -min(vertices(mesh));
            move(mesh,v)
            h = 1/max(vertices(mesh),[],'all');
            scale(mesh,h)
        end
        function vert = vertices(mesh,new_vert)
            if nargin == 1
                vert = get0(mesh,1:mesh.elem_num(1));
            else % nargout = 0 in this case
                set0(mesh,new_vert,1:mesh.elem_num(1));
            end
        end
        function edges = edges(mesh)
            edges = get1(mesh,1:mesh.elem_num(2));
        end
        function dim = manifold(mesh)
            dim = find(mesh.elem_num>0,1,'last')-1;
        end
        function [h,hvect] = meshsize(mesh)
            dim = manifold(mesh);
            hvect = zeros(1,mesh.elem_num(dim+1));
            for i = 1:mesh.elem_num(dim+1)
                elem = RView(mesh,i,dim);
                hvect(i) = max(pdist(vertices(elem)));
            end
            h = max(hvect);
        end
        function view = boundary(mesh)
            dim = manifold(mesh);
            if isempty(dim)
                view = RView(mesh,[],0); % empty mesh
            else
                B_ids = [];
                B_dim = dim-1;
                for i = 1:mesh.elem_num(B_dim+1)
                    if length(get_share(mesh,i,B_dim)) < 2
                       B_ids = [B_ids,i];
                    end
                end
                view = RView(mesh,B_ids,B_dim);
            end
        end
        function plot(mesh,numbers_dim)
            if mesh.dim == 2 && nargin == 1
                plot_edges(mesh)
                return
            end
            
            if nargin ==2
                plot_edges(mesh)
                for i = 1:mesh.elem_num(numbers_dim+1)
                    v = RView(mesh,i,numbers_dim);
                    c = mean(vertices(v),1);
                    if mesh.dim == 3
                        text(c(1),c(2),c(3),num2str(i));
                    else
                        text(c(1),c(2),num2str(i));
                    end
                end
                return
            end
            
            
            hold on
            B = boundary(mesh);
            for i = 1:B.elem_num(3)
                face = RView(B,i,2);
                P = edges2polygon(vertices(face),edges(face));
                fill3(P(:,1),P(:,2),P(:,3),'w','FaceAlpha',1)
%                 fill3(P(:,1),P(:,2),P(:,3),'r','FaceAlpha',0.2)
            end
            axis equal
            axis off
            view(3)
            hold off
            
            
        end

    end

end
