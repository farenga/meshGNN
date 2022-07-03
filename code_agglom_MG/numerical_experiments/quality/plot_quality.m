clear, clc, close all

grids = {'tria','rand_tria','rand_voro','quads'};
methods = {'metis','kmeans','GNN_base','GNN_Res'};
names = {'metis','kmeans','SAGE-Base','SAGE-Res'};

HV = 'V'; % 'H' = horizontal, 'V' = vertical 
format = @(HV,H,V) (HV == 'H')*H + (HV == 'V')*V;

M = length(methods);
G = length(grids);
bins = 10;

switch HV
    case 'H'
        tiledlayout(2,G,'Padding', 'none', 'TileSpacing', 'loose');
    case 'V'
        tiledlayout(G,2,'Padding', 'loose', 'TileSpacing', 'loose');
end
for i = 1:G
    for j = 1:M
        load(['quality_',grids{i},'_',methods{j}])
        %%%
        switch HV
            case 'H'
                nexttile(i)
            case 'V'
                nexttile(2*(i-1)+1)
        end
        hold on
        histogram(UF,bins,'BinLimits',[0,1],'Normalization','probability')
%         ylim([0,1])
        %%%
        switch HV
            case 'H'
                nexttile(i+G)
            case 'V'
                nexttile(2*(i-1)+2)
        end
        hold on
        histogram(CR,bins,'BinLimits',[0,1],'Normalization','probability')
        %%%
%         ylim([0,1])
    end
end

%% labels
font = 12;

grids = {'triangles','random','voronoi','squares'};

for i = 1:G
    switch HV
        case 'H'
            nexttile(i)
            title(grids{i},'fontsize',font)
        case 'V'
            nexttile(2*(i-1)+1)
            ylabel(grids{i},'fontweight','bold','fontsize',font)
    end
    
end

switch HV
    case 'H'
    	nexttile(1)
        ylabel('Uniformity Factor','fontweight','bold','fontsize',font)
        nexttile(G+1)
        ylabel('Circle Ratio','fontweight','bold','fontsize',font)
    case 'V'
        nexttile(1)
        title('Uniformity Factor','fontweight','bold','fontsize',font)
        nexttile(2)
        title('Circle Ratio','fontweight','bold','fontsize',font)
end

f = gcf;
L = legend(names,'Orientation','horizontal',...
    'fontsize',10);

switch HV
    case 'H'
        f.Name = 'quality_horizontal';
        f.WindowState = 'maximized';
        L.Position = [0.4068 0.4605 0.2204 0.0374];
    case 'V'
        f.Name = 'quality_vertical';
        L.Position = [0.24 0.0272 0.56 0.0378];
        f.Position = [489 41.8000 748.8000 740.8000];
end

save_all_figures
