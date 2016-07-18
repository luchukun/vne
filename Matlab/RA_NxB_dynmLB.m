clear all;
close all;

%Fat-tree

N=12;
B=[200,250,300,350,400];
success_opt=[1 1 1 0.986667 0.983333   ];   %success_rate 
success_4=[ 1 1 0.996667 0.993333 0.97 ];   %success_rate
success_2=[1 0.99 0.99 0.776667 0.276667 ];   %success_rate 
success_1=[ 0.97 0.576667 0.06 0 0 ];   %success_rate 
success_4_dyn = [1 1 1 0.98 0.92 ] ;
success_2_dyn =[ 1 0.99 0.95 0.52 0.04  ];
success_1_dyn = [0.98 0.65 0.04 0 0 ]; 

 %1 0.98 0.87 0.56 0.12 

%figure;
subplot(1,2,1)
h1=plot(B,success_opt,B,success_1,B,success_2,B,success_4,B,success_1_dyn,B,success_2_dyn,B,success_4_dyn);
title('(c) Fat-tree, N=12');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'Optimal','LB(K=1)','LB(K=2)','LB(K=4)','dynamicLB(K=1)','dynamicLB(K=2)','dynamicLB(K=4)');
set(gca,'XLim',[200,400]);set(gca,'XGrid','on','YGrid','on')
set(h1,{'Color'},{'r';'b';'g';'k';'b';'g';'k';});
set(h1,{'Marker'},{'*';'O';'s';'+';;'O';'s';'+';},{'MarkerSize'},{10;8;8;10;8;8;10});
set(h1,'LineWidth',2,{'LineStyle'},{'-';'-';'-';'-';':';':';':';});
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%BCube
N=12;
B=[200,250,300,350,400];
success_opt=[  1 1 1 1 1  ];   %success_rate 
success_1=[ 1 0.94 0.45 0 0 ];   %success_rate 
success_2=[ 1 1 1 0.69 0.1  ];   %success_rate 
success_2_dyn =[ 1 1 0.7 0.1 0 ];
success_1_dyn = [1 0.87 0.23 0 0]; 
%success_4=[ 0.999 0.999 0.964 0.533 0.071  ];   %success_rate 
 %1 0.98 0.87 0.56 0.12 

%figure;
subplot(1,2,2);

h1=plot(B,success_opt,B,success_1,B,success_2,B,success_1_dyn,B,success_2_dyn);
title('(d) BCube, N=12');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'Optimal','LB(K=1)','LB(K=2)','dynamicLB(K=1)','dynamicLB(K=2)');
set(gca,'XLim',[200,400]);set(gca,'XGrid','on','YGrid','on')
set(h1,{'Color'},{'r';'b';'g';'b';'g';});
set(h1,{'Marker'},{'*';'O';'s';'O';'s'},{'MarkerSize'},{10;8;8;8;8;});
set(h1,'LineWidth',2,{'LineStyle'},{'-';'-';'-';':';':'});
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 


