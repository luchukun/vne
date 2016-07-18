clear all;
close all;


%Fat-tree
N=[2,4,6,8,9,10];
suc_randrop=[ 1 0.999 0.987 0.858 0.739 0.58 ];   %success_rate backtracking
suc_pert=[ 1 1 0.98 0.875 0.751 0.618];   %success_rate pertuabtion 

figure;
h1=plot(N,suc_randrop,N,suc_pert);
title('Fat-tree');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'ramdom drop','pertubation');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'});
set(h1,'LineWidth',1);

%BCube
N=[2,4,6,8,9,10];
suc_randrop=[1 1 1 0.974 0.882 0.671];   %success_rate backtracking
suc_pert=[ 1 1 1 0.98 0.875 0.682 ];   %success_rate pertuabtion 

figure;
h1=plot(N,suc_randrop,N,suc_pert);
title('BCube');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'ramdom drop','pertubation');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'});
set(h1,'LineWidth',1);

%Vl2
N=[2,4,6,8,9,10];
suc_randrop=[1 1 1 0.97 0.94 0.85 ];   %success_rate backtracking
suc_pert=[  1 1 1 0.98 0.96 0.89  ];   %success_rate pertuabtion 

figure;
h1=plot(N,suc_randrop,N,suc_pert);
title('VL2');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'ramdom drop','pertubation');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'});
set(h1,'LineWidth',1);


