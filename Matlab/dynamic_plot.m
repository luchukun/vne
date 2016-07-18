clear all;
close all;
%VL2
p=[0.8 0.6 0.4 0.2];
success_b=[0.8 0.91 0.96 1 ];   %success_rate backtracking
success_p=[ 0.81 0.86 0.94 0.99];   %success_rate pertuabtion
success_f=[0.77 0.82 0.94 0.97];   %success_rate fisrt fit
figure;
h1=plot(p,success_b,p,success_p,p,success_f);
title('VL2, K=2,N=2~10, B=100~700');

xlabel('Load');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit');
set(gca,'XLim',[0.2,  0.8]);
set(h1,{'Color'},{'r';'b';'g';});
set(h1,{'Marker'},{'*';'O';'s';});
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%Bcube
p=[0.8 0.6 0.4 0.2];
success_b=[0.85 0.94 0.95 0.98 ];   %success_rate backtracking
success_p=[ 0.76 0.84 0.92 0.97];   %success_rate pertuabtion
success_f=[0.74 0.78 0.82 0.89 ];   %success_rate fisrt fit
figure;
h1=plot(p,success_b,p,success_p,p,success_f);
title('Bcube, K=2,N=2~10, B=100~700');

xlabel('Load');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit');
set(gca,'XLim',[0.2,  0.8]);
set(h1,{'Color'},{'r';'b';'g';});
set(h1,{'Marker'},{'*';'O';'s';});
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%FatTree
p=[0.8 0.6 0.4 0.2];
success_b=[0.81 0.91 0.96 1  ];   %success_rate backtracking
success_p=[  0.85 0.9 0.94 1 ];   %success_rate pertuabtion
success_f=[ 0.84 0.9 0.93 1 ];   %success_rate fisrt fit
figure;
h1=plot(p,success_b,p,success_p,p,success_f);
title('Fattree, K=2,N=2~10, B=100~700');

xlabel('Load');
ylabel('Success rate');
legend(h1,'backtracking','pertubation','first-fit');
set(gca,'XLim',[0.2,  0.8]);
set(h1,{'Color'},{'r';'b';'g';});
set(h1,{'Marker'},{'*';'O';'s';});
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 