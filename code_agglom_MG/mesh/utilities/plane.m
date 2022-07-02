function [n,x0,score] = plane(points)
N = size(points,1);
if  N == 2
    x0 = (points(1,:)+points(2,:))/2;
    n = (points(2,:)-points(1,:));
    n = n/norm(n);
elseif N > 2
    x0 = mean(points,1);
%     x0 = points(1,:);
    [coeff,score] = pca(points,'Economy',false);
    n = coeff(:,end)';
else % N < 2
    error('not enought points') 
end

end