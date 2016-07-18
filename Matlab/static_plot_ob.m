clear all;
close all;

%VL2
f=[1,2,4,5];%oversubscription factor
success_opt=[  1 1 0.998 0.972   ];   %success_rate 
success_4=[1 1 0.966 0.936 ];   %success_rate 
success_2=[ 1 0.99 0.956 0.928  ];   %success_rate 
success_1=[ 0.998 0.982 0.809 0.775  ];   %success_rate 

 %1 0.98 0.87 0.56 0.12 

figure;
%subplot(1,2,1);
h1=plot(f,success_opt,f,success_1,f,success_2,f,success_4);
title('VL2,N=8,B=100~700Mbps');
xlabel('Oversubscription ratio');
ylabel('Success rate');
legend(h1,'Optimal','LB(K=1)','LB(K=2)','LB(K=4)');
set(gca,'XLim',[0.5,5.5]);
set(gca,'XTick',[1,2,4,5])
set(gca,'XTickLabel',{'1:1','2:1','4:1','5:1'});set(gca,'XGrid','on','YGrid','on');
set(h1,{'Color'},{'r';'b';'g';'k';});
set(h1,{'Marker'},{'*';'O';'s';'+';},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',2);
%VL2
%N=8;B=100-700;

f=[1,2,4,5];%oversubscription factor


suc_bktrack=[1 1 0.97 0.95  ];   %success_rate backtrack
suc_pert=[1 0.99 0.95 0.93 ];   %success_rate pertuabtion
suc_first=[  1 0.94 0.76 0.66  ];   %success_rate fisrt fit
suc_next=[ 1 0.96 0.77 0.71 ];   %success_rate fisrt fit
suc_best=[ 0.99 0.92 0.75 0.66  ];   %success_rate fisrt fit 
figure;
%subplot(1,2,2);
h1=plot(f,suc_bktrack,f,suc_pert,f,suc_first,f,suc_next,f,suc_best);
title('VL2,N=8,B=100~700Mbps');
xlabel('Oversubscription ratio');
ylabel('Success rate');

legend(h1,'backtrack','pertubation','first-fit','next-fit','greedy','Location','SouthWest');
set(gca,'XLim',[0.5,5.5]);
set(gca,'XTick',[1,2,4,5])
set(gca,'XTickLabel',{'1:1','2:1','4:1','5:1'});set(gca,'XGrid','on','YGrid','on')
set(h1,{'Color'},{'r';'b';'g';'m';'c'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'},{'MarkerSize'},{10;8;8;10;10});
set(h1,'LineWidth',2);

