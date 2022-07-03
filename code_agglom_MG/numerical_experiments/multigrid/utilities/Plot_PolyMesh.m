function [error] = Plot_PolyMesh(Tau,id_plot,color,neighbour,line,color_boundary,linewidth)

if(nargin<7)
    linewidth=2;
end

error = 1;
figure(id_plot);
if color == '0'
    mark_boundary = [0.76          0.97          0.04];
    mark = [0.76          0.97          0.04];
    coord = Tau.coords_element{1};
    plot([coord(:,1); coord(1,1)], [coord(:,2); coord(1,2)], 'Color', mark, 'Linewidth',linewidth);
    hold on; axis equal;
    for ie = 2:Tau.ne
        coord = Tau.coords_element{ie};
        coord = [coord; coord(1,:)];
        for iedg = 1:neighbour.nedges(ie)
            if neighbour.neigh{ie}(iedg) == -1
                plot([coord(iedg,1) coord(iedg+1,1)], [coord(iedg,2) coord(iedg+1,2)], 'Color', mark_boundary, 'Linewidth',linewidth);
            else
                plot([coord(iedg,1) coord(iedg+1,1)], [coord(iedg,2) coord(iedg+1,2)], 'Color', mark, 'Linewidth',linewidth);
            end
        end
    end
else
    mark_boundary = [line,color_boundary];
    mark = [line,color];
    coord = Tau.coords_element{1};
    %plot([coord(1,1)], [coord(1,2)], mark, 'Linewidth',2);
    hold on; axis equal;
    for ie = 1:Tau.ne
        coord = Tau.coords_element{ie};
        coord = [coord; coord(1,:)];
        for iedg = 1:neighbour.nedges(ie)
            if neighbour.neigh{ie}(iedg) == -1
                plot([coord(iedg,1) coord(iedg+1,1)], [coord(iedg,2) coord(iedg+1,2)], mark_boundary, 'Linewidth',linewidth);
            else
                plot([coord(iedg,1) coord(iedg+1,1)], [coord(iedg,2) coord(iedg+1,2)], mark, 'Linewidth',linewidth);
            end
        end
    end
end
axis off;


%figure(id_plot+1);
%mark_boundary = ['-',color_boundary];
%mark = ['-',color];
%coord = Tau.coords_element{1};
%plot([coord(:,1); coord(1,1)], [coord(:,2); coord(1,2)], mark, 'Linewidth',2);


% BBox = Tau.BBox(1,:);
% x1B = BBox(1); x2B = BBox(2);
% y1B = BBox(3); y2B = BBox(4);
% [X,Y] = meshgrid(x1B:(x2B-x1B)/10:x2B, y1B:(y2B-y1B)/10:y2B);
% func = ones(size(X,1),size(X,2));
% surf(X,Y,func);
% hold on; axis equal;
% for ie = 2:Tau.ne
%     BBox = Tau.BBox(ie,:);
%     x1B = BBox(1); x2B = BBox(2);
%     y1B = BBox(3); y2B = BBox(4);
%     [X,Y] = meshgrid(x1B:(x2B-x1B)/10:x2B, y1B:(y2B-y1B)/10:y2B);
%     if mod(ie,2)==0
%         func = zeros(size(X,1),size(X,2));
%     else
%         func = ones(size(X,1),size(X,2));
%     end
%     surf(X,Y,func);
% end


error = 0;

return