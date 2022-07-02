function plot_edges(mesh,clust)
% plot mesh wire-frame

hold on

if nargin == 2
    labels = unique(clust);
    n_clust = length(labels);
    my_legend = cell(1,n_clust);
    switch mesh.dim
        case 1
            for i = 1:n_clust
                my_legend{i} = ['cluster ',num2str(i)];
                x = get0(mesh,find(clust == labels(i)));
                plot(x,zeros(size(x)),'.','MarkerSize',16)
            end
        case 2
            for i = 1:n_clust
                my_legend{i} = ['cluster ',num2str(i)];
                vert = get0(mesh,find(clust == labels(i)));
                plot(vert(:,1),vert(:,2),'.','MarkerSize',16)
            end
        otherwise
            for i = 1:n_clust
                my_legend{i} = ['cluster ',num2str(i)];
                vert = get0(mesh,find(clust == labels(i)));
                plot3(vert(:,1),vert(:,2),vert(:,3),'.','MarkerSize',16)
            end
    end

    legend(my_legend,'AutoUpdate','off')
end


axis equal
axis off

switch mesh.dim
    case 1
        for i = 1:mesh.elem_num(2)
            edge = get1(mesh,i);
            x = get0(mesh,edge);
            plot(x,zeros(size(x)),'-k')
        end
    case 2
        for i = 1:mesh.elem_num(2)
            edge = get1(mesh,i);
            vert = get0(mesh,edge);
            plot(vert(:,1),vert(:,2),'-k')
        end
    otherwise
        for i = 1:mesh.elem_num(2)
            edge = get1(mesh,i);
            vert = get0(mesh,edge);
            plot3(vert(:,1),vert(:,2),vert(:,3),'-k')
        end
        view(3)
end

hold off