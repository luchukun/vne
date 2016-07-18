clear all;
close all;
%VL2
N=[2,4,6,8,9,10];


suc_bktrack=[1 1 1 0.975 0.965 0.895 ];   %success_rate backtracking
suc_pert=[1 1 0.992 0.968 0.92 0.839 ];   %success_rate pertuabtion
suc_first=[0.999 0.987 0.817 0.531 0.424 0.321 ];   %success_rate fisrt fit
suc_next=[0.999 0.985 0.845 0.58 0.477 0.38 ];   %success_rate fisrt fit
suc_best=[0.998 0.984 0.824 0.595 0.489 0.413 ];   %success_rate fisrt fit 
figure;
h1=plot(N,suc_bktrack,N,suc_pert,N,suc_first,N,suc_next,N,suc_best);
title('VL2');
xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'},{'MarkerSize'},{10;8;8;10;10});
set(h1,'LineWidth',1.5);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%Fat-tree
N=[2,4,6,8,9,10];
suc_bktrack=[1 1 1 0.88 0.83 0.65 ];   %success_rate backtracking
suc_pert=[ 1 1 0.98 0.875 0.751 0.618];   %success_rate pertuabtion 
suc_first=[ 0.999 0.985 0.792 0.455 0.308 0.182 ];   %success_rate fisrt fit
suc_next=[  0.999 0.985 0.811 0.468 0.348 0.238  ];   %success_rate fisrt fit
suc_best=[ 0.998 0.983 0.803 0.512 0.383 0.28 ];   %success_rate fisrt fit 
figure;
h1=plot(N,suc_bktrack,N,suc_pert,N,suc_first,N,suc_next,N,suc_best);
title('Fat-tree');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'},{'MarkerSize'},{10;8;8;10;10});
set(h1,'LineWidth',1.5);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 


%Bcube
N=[2,4,6,8,9,10];
suc_bktrack=[  1 1 1 0.995 0.97 0.83  ];   %success_rate backtracking
suc_pert=[1 1 1 0.995 0.895 0.69 ];   %success_rate pertuabtion
suc_first=[ 1 0.995 0.885 0.54 0.315 0.19 ];   %success_rate fisrt fit
suc_next=[1 1 0.905 0.55 0.34 0.205 ];   %success_rate fisrt fit
suc_best=[ 1 0.995 0.935 0.72 0.475 0.285 ];   %success_rate fisrt fit 
figure;
h1=plot(N,suc_bktrack,N,suc_pert,N,suc_first,N,suc_next,N,suc_best);
title('BCube');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'},{'MarkerSize'},{10;8;8;10;10});
set(h1,'LineWidth',1.5);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%Tree
N=[2,4,6,8,9,10];
suc_bktrack=[  1 1 0.985 0.924 0.851 0.73 ];   %success_rate backtracking
suc_pert=[ 1 0.999 0.983 0.902 0.803 0.697  ];   %success_rate pertuabtion
suc_first=[  1 0.992 0.868 0.654 0.467 0.36 ];   %success_rate fisrt fit
suc_ar=[  1 1 0.882 0.526 0.278 0.164 ];   %success_rate allocation range
%20servers and 1Gbps
figure;
h1=plot(N,suc_bktrack,N,suc_pert,N,suc_first,N,suc_ar);
title('Tree');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','allocation-range','Location','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'d'},{'MarkerSize'},{10;8;8;8});
set(h1,'LineWidth',1.5);



