close all;
linewidth = 3;
markersize = 5;
marker = 'ox^*sp'; marker = [marker, marker];
colors = 'rbkmcg'; colors = [colors, colors];
id = 0;

plot_nested = 1;
plot_non_nested = 1;

x = 1:9;
axis_vect_cond = [x(1)-0.25 10 20 6*10^3]; yticks_vec_cond = [30 100 500 1000 5000];
axis_vect_iter = [x(1)-0.25 10 20 8*10^2]; yticks_vec_iter = [30 80 300 700];

C_p1_c = 30;
C_p2 = 80;

C_p05 = 40;
C_p1_i = 90;

%% CONDITION NUMBER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT CONDITION NUMBERS   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = 1:9;

if plot_non_nested
    
    %NON NESTED 262vs17
    % id = id + 1;
    % line = [colors(id),marker(id),'--'];
    % cond = [158.4495 412.0293 875.4619 1555.308 2278.6017 3170.8677 4142.9754 5163.978 6992.0863];
    % loglog(x,cond,line,'LineWidth',3,'MarkerSize',15); hold on;
    % legendInfo{id} = ['Non-nest. poly. N_h=216 N_H=17'];
    
    %NON NESTED Tau_h=262 vs Tau_H=18 where Tau_H=18 is Voronoi
    id = id + 1;
    line = [colors(id),marker(id),'-.'];
    cond = [60.5186 209.6212 488.3035 813.926 1317.8976 1985.7064 2592.3762 3238.7525 3905.5234];
    loglog(x,cond,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Non-nest. poly. N_h=216'];
    
    %NON NESTED Tau_h=516 vs Tau_H=32
    % id = id + 1;
    % line = [colors(id),marker(id),'-.'];
    % cond = [95.961 477.4693 1012.1561 1726.5552 2569.7184 3445.2898 4426.0296 5622.9991 6859.0039];
    % loglog(x,cond,line,'LineWidth',3,'MarkerSize',15); hold on;
    % legendInfo{id} = ['Non-nest. poly. N_h=516 N_H=32'];
    
    %NON NESTED Tau_h=262 vs Tau_H=4, where Tau_H=4 is Voronoi
    % id = 1;
    % line = [colors(id),marker(id),'-.'];
    % cond = [128.1817 350.9397 670.6553 1107.3387 1593.6165 2146.4288 2678.2044 3284.5908 3337.483];
    % loglog(x,cond,line,'LineWidth',3,'MarkerSize',15); hold on;
    % legendInfo{id} = ['Non-nest. poly. N_h=216 N_H=4'];
    
    %NON NESTED Tau_h=516 vs Tau_H=16, where Tau_H=16 is Voronoi
    %id = id + 1;
    %line = [colors(id),marker(id),'-.'];
    %cond = [87.3992 312.8201 684.5727 1136.5796 1571.4078 2101.9069 2571.4514 3083.6588 3869.5449];
    %loglog(x,cond,line,'LineWidth',3,'MarkerSize',15); hold on;
    %legendInfo{id} = ['Non-nest. poly. N_h=516'];
    
    %NON NESTED Tau_h=516 vs Tau_H=64, where Tau_H=64 is Voronoi
    id = id + 1;
    line = [colors(id),marker(id),'-.'];
    cond = [34.3035 178.7451 418.7812 682.6494 1060.6766 1534.0321 2101.8662 2776.7924 3597.7841];
    loglog(x,cond,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Non-nest. poly. N_h=516'];
end


if plot_nested
    %NESTED Tau_h=262 vs Tau_H=17
    id = id + 1;
    line = [colors(id),marker(id),'-'];
    cond = [62.641 160.0706 241.8657 325.4846 379.2382 456.2464 551.1343 603.7709 661.4906];
    loglog(x,cond,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Nested poly. N_h=216'];
    
    %NESTED Tau_h=516 vs Tau_H=34
    id = id + 1;
    line = [colors(id),marker(id),'-'];
    cond = [64.9171 179.6452 288.2962 385.0758 483.7439 533.9126 625.9311 624.4185 675.1565];
    loglog(x,cond,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Nested poly. N_h=516'];
    
    %NESTED QUADRILATERAL Tau_h=256 vs Tau_H=64
    id = id + 1;
    line = [colors(id),marker(id),'-'];
    cond = [43.4183 115.4801 194.2126 252.5511 336.7532 411.7016 543.2502 496.9999 631.2337];
    loglog(x,cond,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Nested quad. N_h=256'];
    
    %NESTED QUADRILATERAL Tau_h=1024 vs Tau_H=64
    id = id + 1;
    line = [colors(id),marker(id),'-'];
    cond = [43.6818 115.4801 194.2126 252.5511 336.7532 411.7016 543.2502 501.3266 633.7908];
    loglog(x,cond,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Nested quad. N_h=1024'];
end

id = id + 1;
loglog(x,C_p1_c*x,'k--','LineWidth',2.5);
legendInfo{id} = ['O(p)'];

id = id + 1;
loglog(x,C_p2*x.^2,'k-.','LineWidth',2.5);
legendInfo{id} = ['O(p^2)'];

grid on;
LEG = legend(legendInfo,'Location','SouthEast'); clear legendInfo;
XLA = xlabel('p');
YLA = ylabel('K(P_{ad})');
TIT = title('Condition number');
set(gca,'fontsize',10,'fontweight','bold');
axis(axis_vect_cond);
yticks(yticks_vec_cond);

%% ITERATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   PLOT ITERATION COUNTS    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
id = 0;
x = 1:9;

if plot_non_nested
    %NON NESTED Tau_h=262 vs Tau_H=17
    % id = id + 1;
    % line = [colors(id),marker(id),'-.'];
    % iter = [97 174 236 310 373 431 484 540 701];
    % loglog(x,iter,line,'LineWidth',3,'MarkerSize',15); hold on;
    % legendInfo{id} = ['Non-nest. poly. N_h=262'];
    
    %NON NESTED Tau_h=262 vs Tau_H=18 where Tau_H=18 is Voronoi
    id = id + 1;
    line = [colors(id),marker(id),'-.'];
    iter = [67 122 170 213 264 314 360 410 472];
    loglog(x,iter,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Non-nest. poly. N_h=216'];
    
    %NON NESTED Tau_h=516 vs Tau_H=32
    % id = id + 1;
    % line = [colors(id),marker(id),'-.'];
    % iter = [107 198 292 372 439 515 590 669 781];
    % loglog(x,iter,line,'LineWidth',3,'MarkerSize',15); hold on;
    % legendInfo{id} = ['Non-nest. poly. N_h=516'];
    
    
    %NON NESTED Tau_h=262 vs Tau_H=4, where Tau_H=4 is Voronoi
    % id = id + 1;
    % line = [colors(id),marker(id),'-.'];
    % iter = [98 181 236 283 327 367 406 457 602];
    % loglog(x,iter,line,'LineWidth',3,'MarkerSize',15); hold on;
    % legendInfo{id} = ['Non-nest. poly. N_h=216'];
    
    %NON NESTED Tau_h=516 vs Tau_H=16, where Tau_H=16 is Voronoi
    %id = id + 1;
    %line = [colors(id),marker(id),'--.'];
    %iter = [88 151 202 244 277 313 354 391 489];
    %loglog(x,iter,line,'LineWidth',3,'MarkerSize',15); hold on;
    %legendInfo{id} = ['Non-nest. poly. N_h=516'];
    
    %NON NESTED Tau_h=516 vs Tau_H=64, where Tau_H=64 is Voronoi
    id = id + 1;
    line = [colors(id),marker(id),'-.'];
    iter = [54 112 172 218 267 320 375 431 494];
    loglog(x,iter,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Non-nest. poly. N_h=516'];
end

if plot_nested
    % NESTED Tau_h=262 vs Tau_H=17
    id = id + 1;
    line = [colors(id),marker(id),'-'];
    iter = [76 120 140 157 167 179 190 197 230];
    loglog(x,iter,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Nested poly. N_h=262'];
    
    % NESTED Tau_h=516 vs Tau_H=34
    id = id + 1;
    line = [colors(id),marker(id),'-'];
    iter = [78 117 138 152 162 169 181 196 232];
    loglog(x,iter,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Nested poly. N_h=516'];
    
    % NESTED QUADRILATERAL Tau_h=256 vs Tau_H=16
    id = id + 1;
    line = [colors(id),marker(id),'-'];
    iter = [55 94 121 136 151 164 175 173 185];
    loglog(x,iter,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Nested quad. N_h=256'];
    
    % NESTED QUADRILATERAL Tau_h=1024 vs Tau_H=64
    id = id + 1;
    line = [colors(id),marker(id),'-'];
    cond = [61 101 121 136 150 161 176 175 188];
    loglog(x,cond,line,'LineWidth',3,'MarkerSize',15); hold on;
    legendInfo{id} = ['Nested quad. N_h=1024'];
end

id = id + 1;
loglog(x,C_p05*sqrt(x),'k--','LineWidth',2.5);
legendInfo{id} = ['O(p^{0.5})'];

id = id + 1;
loglog(x,C_p1_i*x,'k-.','LineWidth',2.5);
legendInfo{id} = ['O(p)'];

grid on;
LEG = legend(legendInfo,'Location','SouthEast'); clear legendInfo;
XLA = xlabel('p');
YLA = ylabel('iter(P_{ad})');
TIT = title('Iteration counts');
set(gca,'fontsize',10,'fontweight','bold');
axis(axis_vect_iter);
yticks(yticks_vec_iter);