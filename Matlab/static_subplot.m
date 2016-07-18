clear all;
close all;
%Fat-tree
N=[2,4,6,8,9,10];
suc_bktrack=[1 1 1 0.88 0.83 0.65 ];   %success_rate backtracking
suc_pert=[ 1 1 0.98 0.875 0.751 0.618];   %success_rate pertuabtion 
suc_first=[ 0.999 0.985 0.792 0.455 0.308 0.182 ];   %success_rate fisrt fit
suc_next=[  0.999 0.985 0.811 0.468 0.348 0.238  ];   %success_rate fisrt fit
suc_best=[ 0.998 0.983 0.803 0.512 0.383 0.28 ];   %success_rate fisrt fit 
subplot(3,1,1);

h1=plot(N,suc_bktrack,N,suc_pert,N,suc_first,N,suc_next,N,suc_best);
title('(a) Fat-tree');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'});
set(h1,'LineWidth',1);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 
%VL2
N=[2,4,6,8,9,10];


suc_bktrack=[1 1 1 0.975 0.965 0.895 ];   %success_rate backtracking
suc_pert=[1 1 0.992 0.968 0.92 0.839 ];   %success_rate pertuabtion
suc_first=[0.999 0.987 0.817 0.531 0.424 0.321 ];   %success_rate fisrt fit
suc_next=[0.999 0.985 0.845 0.58 0.477 0.38 ];   %success_rate fisrt fit
suc_best=[0.998 0.984 0.824 0.595 0.489 0.413 ];   %success_rate fisrt fit 
subplot(3,1,2);
h1=plot(N,suc_bktrack,N,suc_pert,N,suc_first,N,suc_next,N,suc_best);
title('(b) VL2');
xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'});
set(h1,'LineWidth',1);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%Bcube
N=[2,4,6,8,9,10];
suc_bktrack=[  1 1 1 0.995 0.97 0.83  ];   %success_rate backtracking
suc_pert=[1 1 1 0.995 0.895 0.69 ];   %success_rate pertuabtion
suc_first=[ 1 0.995 0.885 0.54 0.315 0.19 ];   %success_rate fisrt fit
suc_next=[1 1 0.905 0.55 0.34 0.205 ];   %success_rate fisrt fit
suc_best=[ 1 0.995 0.935 0.72 0.475 0.285 ];   %success_rate fisrt fit 
subplot(3,1,3);

h1=plot(N,suc_bktrack,N,suc_pert,N,suc_first,N,suc_next,N,suc_best);
title('(c) BCube');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'});
set(h1,'LineWidth',1);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%Tree
N=[2,4,6,8,9,10];
suc_bktrack=[ 1 0.998 0.973 0.896 0.811 0.721 ];   %success_rate backtracking
suc_pert=[ 1 0.998 0.978 0.862 0.779 0.654 ];   %success_rate pertuabtion
suc_first=[ 1 0.984 0.861 0.618 0.499 0.376 ];   %success_rate fisrt fit
suc_next=[ 0.999 0.964 0.794 0.536 0.412 0.323  ];   %success_rate fisrt fit
suc_best=[ 1 0.982 0.824 0.561 0.454 0.338  ];   %success_rate fisrt fit 
suc_ar=[  1 1 0.845 0.487 0.297 0.164  ];   %success_rate allocation range

figure;
h1=plot(N,suc_bktrack,N,suc_pert,N,suc_first,N,suc_next,N,suc_best,N,suc_ar);
title('Tree');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','allocation-range','Location','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'m';'c';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x';'d'});
set(h1,'LineWidth',1);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 


