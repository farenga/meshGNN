function plot_grids(grids,boom)
if nargin < 2
    boom = false;
end
[r,c] = size(grids);
tiledlayout(r,c, 'Padding', 'tight', 'TileSpacing', 'none');

for i = 1:r
    for j = 1:c
        id = (i-1)*c+j;
        nexttile(id)
        if boom
            explode(grids{i,j})
        else
            plot(grids{i,j})
        end
    end
end