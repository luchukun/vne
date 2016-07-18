clear all;
close all;


%Fat-tree
N=[2,4,6,8,9,10];
suc_bktrack=[1 1 1 0.88 0.83 0.65 ];   %success_rate backtracking
suc_pert=[ 1 1 0.98 0.875 0.751 0.618];   %success_rate pertuabtion 
suc_first=[ 0.999 0.985 0.792 0.455 0.308 0.182 ];   %success_rate fisrt fit
suc_next=[  0.999 0.985 0.811 0.468 0.348 0.238  ];   %success_rate fisrt fit
suc_best=[ 0.998 0.983 0.803 0.512 0.383 0.28 ];   %success_rate fisrt fit 
figure;
subplot(1,4,1);
h1=plot(N,suc_bktrack,N,suc_pert,N,suc_first,N,suc_next,N,suc_best);
title('(a) Fat-tree,B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'},{'MarkerSize'},{10;8;8;10;10});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 
%Fat-tree
N=8;
B=[200 250 300 350 400] ;
suc_bktrack=[ 1 1 1 0.99 0.92  ];   %success_rate backtracking
suc_pert=[ 1 1 1 0.97 0.92   ];   %success_rate pertuabtion
suc_first=[ 0.93 0.92 0.87 0.66 0.37  ];   %success_rate fisrt fit
suc_next=[ 0.96 0.93 0.9 0.69 0.46  ];   %success_rate fisrt fit
suc_best=[  0.95 0.89 0.81 0.65 0.45  ];   %success_rate fisrt fit 
%figure;
subplot(1,4,2);
h1=plot(B,suc_bktrack,B,suc_pert,B,suc_first,B,suc_next,B,suc_best);
title('(b) Fat-tree,N=8');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'},{'MarkerSize'},{10;8;8;10;10});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%Bcube
N=[2,4,6,8,9,10];
suc_bktrack=[  1 1 1 0.995 0.97 0.83  ];   %success_rate backtracking
suc_pert=[1 1 1 0.995 0.895 0.69 ];   %success_rate pertuabtion
suc_first=[ 1 0.995 0.885 0.54 0.315 0.19 ];   %success_rate fisrt fit
suc_next=[1 1 0.905 0.55 0.34 0.205 ];   %success_rate fisrt fit
suc_best=[ 1 0.995 0.935 0.72 0.475 0.285 ];   %success_rate fisrt fit 
%figure;
subplot(1,4,3);
h1=plot(N,suc_bktrack,N,suc_pert,N,suc_first,N,suc_next,N,suc_best);
title('(c) BCube,B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'},{'MarkerSize'},{10;8;8;10;10});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 
%Bcube
N=8;
B=[200 250 300 350 400] ;
suc_bktrack=[  1 1 1 1 1 ];   %success_rate backtracking
suc_pert=[ 1 1 1 1 0.98  ];   %success_rate pertuabtion
suc_first=[ 0.96 0.95 0.81 0.59 0.44 ];   %success_rate fisrt fit
suc_next=[0.97 0.94 0.82 0.61 0.46 ];   %success_rate fisrt fit
suc_best=[ 1 0.99 0.9 0.79 0.64 ];   %success_rate fisrt fit 
%figure;
subplot(1,4,4);
h1=plot(B,suc_bktrack,B,suc_pert,B,suc_first,B,suc_next,B,suc_best);
title('(d) BCube,N=8');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'},{'MarkerSize'},{10;8;8;10;10});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 
%Tree
N=[2,4,6,8,9,10];
suc_bktrack=[  1 1 0.985 0.924 0.851 0.73 ];   %success_rate backtracking
suc_pert=[ 1 0.999 0.983 0.902 0.803 0.697  ];   %success_rate pertuabtion
suc_first=[  1 0.992 0.868 0.654 0.467 0.36 ];   %success_rate fisrt fit
suc_ar=[  1 1 0.882 0.526 0.278 0.164 ];   %success_rate allocation range


%Tree
N=[2,4,6,8,9,10];
suc_bktrack=[ 1 1 0.97 0.76 0.62 0.46  ];   %success_rate backtracking
suc_pert=[ 1 0.998 0.945 0.741 0.595 0.39  ];   %success_rate pertuabtion
suc_first=[  0.999 0.969 0.753 0.436 0.312 0.208 ];   %success_rate fisrt fit
suc_ar=[  1 1 0.686 0.234 0.114 0.047 ];   %success_rate allocation range
%20servers and 1Gbps
figure;
subplot(1,2,2);
h1=plot(N,suc_bktrack,N,suc_pert,N,suc_first,N,suc_ar);
title('(b) Tree, B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','allocation-range','Location','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'d'},{'MarkerSize'},{10;8;8;8});
set(h1,'LineWidth',2);

%Tree

B=[200 250 300 350 400] ;
suc_bktrack=[ 0.99 0.99 0.99 0.93 0.79  ];   %success_rate backtracking
suc_pert=[ 0.999 0.998 0.975 0.867 0.647 ];   %success_rate pertuabtion
suc_first=[ 0.918 0.861 0.72 0.508 0.316  ];   %success_rate fisrt fit
suc_ar=[ 0.999 0.985 0.892 0.742 0.586  ];   %success_rate allocation range
%20servers and 1Gbps
%figure;
subplot(1,2,1);
h1=plot(B,suc_bktrack,B,suc_pert,B,suc_first,B,suc_ar);
title('(a) Tree, N=8');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit','allocation-range','Location','SouthWest');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'d'},{'MarkerSize'},{10;8;8;8});
set(h1,'LineWidth',2);

