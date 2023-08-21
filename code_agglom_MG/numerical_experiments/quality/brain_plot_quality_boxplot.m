clear, clc, close all

methods = {'metis','kmeans','GNN'};
names = {'metis','k-means','GNN'};

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


    UF_tot = [];
    CR_tot = [];
    g = [];
    for j = 1:M
        load(['quality_brain_',methods{j}])
        UF_tot = [UF_tot;UF'];
        CR_tot = [CR_tot;CR'];
        g = [g;repmat(names(j),length(UF),1)];
    end
    %%%
        switch HV
        case 'H'
            nexttile(1)
        case 'V'
            nexttile(2*(1-1)+1)
    end
        hold on
        boxplot(UF_tot,g)
        switch HV
        case 'H'
            nexttile(1+1)
        case 'V'
            nexttile(2*(1-1)+2)
    end
        hold on
        boxplot(CR_tot,g)


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


switch HV
    case 'H'
        f.Name = 'quality_horizontal_brain_box';
        f.WindowState = 'maximized';
    case 'V'
        f.Name = 'quality_vertical_brain_box';
        f.Position = [489 41.8000 748.8000 740.8000];
end

save_all_figures
