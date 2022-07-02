clear, clc, close all

methods = {'metis','kmeans','GNN_base','GNN_Res'};
names = {'metis','kmeans','SAGE-Base','SAGE-Res'};

HV = 'V'; % 'H' = horizontal, 'V' = vertical 
format = @(HV,H,V) (HV == 'H')*H + (HV == 'V')*V;

M = length(methods);
bins = 10;

switch HV
    case 'H'
        tiledlayout(2,1,'Padding', 'none', 'TileSpacing', 'loose');
    case 'V'
        tiledlayout(1,2,'Padding', 'loose', 'TileSpacing', 'loose');
end

for j = 1:M
    load(['quality_brain_',methods{j}])
    %%%
    switch HV
        case 'H'
            nexttile(1)
        case 'V'
            nexttile(2*(1-1)+1)
    end
    hold on
    histogram(UF,bins,'BinLimits',[0,1],'Normalization','probability')
%         ylim([0,1])
    %%%
    switch HV
        case 'H'
            nexttile(1+1)
        case 'V'
            nexttile(2*(1-1)+2)
    end
    hold on
    histogram(CR,bins,'BinLimits',[0,1],'Normalization','probability')
    %%%
%         ylim([0,1])
end

%% labels
font = 12;

switch HV
    case 'H'
        nexttile(1)
        title(grids{1},'fontsize',font)
    case 'V'
        nexttile(2*(1-1)+1)
        ylabel('Brain','fontweight','bold','fontsize',font)
end

switch HV
    case 'H'
    	nexttile(1)
        ylabel('Uniformity Factor','fontweight','bold','fontsize',font)
        nexttile(1+1)
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
        f.Name = 'quality_horizontal_brain';
        f.WindowState = 'maximized';
        L.Position = [0.4068 0.4605 0.2204 0.0374];
    case 'V'
        f.Name = 'quality_vertical_brain';
        L.Position = [0.24 0.02 0.56 0.035];
        f.Position = [489 41.8000 748.8000 740.8000];
end

save_all_figures
