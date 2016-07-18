clear all;
close all;

%Fat-tree
N=8;
B=[200 250 300 350 400] ;
suc_bktrack=[ 1 1 1 0.99 0.92  ];   %success_rate backtracking
suc_pert=[ 1 1 1 0.97 0.91 ];   %success_rate pertuabtion
suc_first=[  0.925 0.875 0.8 0.605 0.37 ];   %success_rate fisrt fit
suc_next=[0.94 0.885 0.805 0.625 0.455  ];   %success_rate fisrt fit
suc_best=[  0.925 0.875 0.725 0.615 0.425  ];   %success_rate fisrt fit 
figure;
h1=plot(B,suc_bktrack,B,suc_pert,B,suc_first,B,suc_next,B,suc_best);
title('Fat-tree,N=8');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'},{'MarkerSize'},{10;8;8;10;10});
set(h1,'LineWidth',1.5);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 


%Bcube
N=8;
B=[200 250 300 350 400] ;
suc_bktrack=[  1 1 1 1 1 ];   %success_rate backtracking
suc_pert=[ 1 1 1 1 0.98  ];   %success_rate pertuabtion
suc_first=[ 0.96 0.95 0.81 0.59 0.44 ];   %success_rate fisrt fit
suc_next=[0.97 0.94 0.82 0.61 0.46 ];   %success_rate fisrt fit
suc_best=[ 1 0.99 0.9 0.79 0.64 ];   %success_rate fisrt fit 
figure;
h1=plot(B,suc_bktrack,B,suc_pert,B,suc_first,B,suc_next,B,suc_best);
title('BCube,N=8');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[200,400]);
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



