clear all;
close all;


%Fat-tree
N=[2,4,6,8,9,10];%B=100-700Mbps;
success_opt=[  1 1 1 0.996667 0.996667 1   ];   %success_rate 
success_4=[ 1 1 1 0.996667 0.993333 0.986667  ];   %success_rate
success_2=[  1 1 0.993333 0.963333 0.9 0.773333  ];   %success_rate
success_1=[1 1 0.836667 0.4 0.22 0.0733333 ];   %success_rate 
 
 
 %1 0.98 0.87 0.56 0.12 

figure;
subplot(1,4,1)
h1=plot(N,success_opt,N,success_1,N,success_2,N,success_4);
title('(a) Fat-tree, B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'Optimal','LB(K=1)','LB(K=2)','LB(K=4)');
set(gca,'XLim',[2,10]);set(gca,'XGrid','on','YGrid','on')
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 




%BCube
N=[2,4,6,8,9,10];%B=100-700;
success_opt=[ 1 1 1 1 1 1 ];   %success_rate 
success_1=[1 1 1 0.64 0.47 0.25 ];   %success_rate 
success_2=[ 1 1 1 0.98 0.9 0.75 ];   %success_rate 
%success_4=[  1 1 1 0.99 0.93 0.743 ];   %success_rate 
 %1 0.98 0.87 0.56 0.12 

%figure;
subplot(1,4,2);
h1=plot(N,success_opt,N,success_1,N,success_2);
title('(b) BCube, B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'Optimal','LB(K=1)','LB(K=2)');
set(gca,'XLim',[2,10]);set(gca,'XGrid','on','YGrid','on')
set(h1,{'Color'},{'r';'b';'g'});
set(h1,{'Marker'},{'*';'O';'s';},{'MarkerSize'},{10;8;8;});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 
%Fat-tree

N=12;
B=[200,250,300,350,400];
success_opt=[1 1 1 0.986667 0.983333   ];   %success_rate 
success_4=[ 1 1 0.996667 0.993333 0.97 ];   %success_rate
success_2=[1 0.99 0.99 0.776667 0.276667 ];   %success_rate 
success_1=[ 0.97 0.576667 0.06 0 0 ];   %success_rate 

 %1 0.98 0.87 0.56 0.12 

%figure;
subplot(1,4,3)
h1=plot(B,success_opt,B,success_1,B,success_2,B,success_4);
title('(c) Fat-tree, N=12');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'Optimal','LB(K=1)','LB(K=2)','LB(K=4)');
set(gca,'XLim',[200,400]);set(gca,'XGrid','on','YGrid','on')
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%BCube
N=12;
B=[200,250,300,350,400];
success_opt=[  1 1 1 1 1  ];   %success_rate 
success_1=[ 1 0.94 0.45 0 0 ];   %success_rate 
success_2=[ 1 1 1 0.69 0.1  ];   %success_rate 
%success_4=[ 0.999 0.999 0.964 0.533 0.071  ];   %success_rate 
 %1 0.98 0.87 0.56 0.12 

%figure;
subplot(1,4,4);

h1=plot(B,success_opt,B,success_1,B,success_2);
title('(d) BCube, N=12');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'Optimal','LB(K=1)','LB(K=2)');
set(gca,'XLim',[200,400]);set(gca,'XGrid','on','YGrid','on')
set(h1,{'Color'},{'r';'b';'g'});
set(h1,{'Marker'},{'*';'O';'s';},{'MarkerSize'},{10;8;8;});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 


