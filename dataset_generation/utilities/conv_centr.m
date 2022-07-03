function C = conv_centr(points,egs)

if size(points,2) == 2
    points = edges2polygon(points,egs);
    x = points(:,1); y = points(:,2);
    % fare il sorting dei vertici!!
    pgon = polyshape(x,y,'KeepCollinearPoints',true);
    [x_c,y_c] = centroid(pgon);
    C = [x_c,y_c];
else
    points = points(unique(convhull(points)),:);
    T = delaunay(points);
    N_tri = size(T,1);
    W = zeros(N_tri,1);
    C = 0;
    for m = 1:N_tri
        P = points(T(m,:),:);
        W(m) = 1/6*abs(dot(cross(P(2,:)-P(1,:),P(3,:)-P(1,:)),P(4,:)-P(1,:)));
        C = C + W(m) * mean(P);
    end
    C = C./sum(W);
end