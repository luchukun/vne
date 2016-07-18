clear all;
close all;


%Fat-tree
N=[2,4,6,8,9,10];%B=100-700Mbps;
success_opt=[  1 1 1 1 0.995 0.99  ];   %success_rate 
success_1=[1 1 0.84 0.42 0.23 0.085 ];   %success_rate 
success_2=[  1 1 0.995 0.96 0.925 0.75  ];   %success_rate 
success_4=[ 1 1 1 0.99 0.975 0.905  ];   %success_rate 
 %1 0.98 0.87 0.56 0.12 

figure;
subplot(1,4,1)
h1=plot(N,success_opt,N,success_1,N,success_2,N,success_4);
title('(a) Fat-tree, B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'Optimal','LB(K=1)','LB(K=2)','LB(K=4)','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%Fat-tree

N=12;
B=[200,250,300,350,400];
success_opt=[1 0.998 0.999 0.993 0.983 ];   %success_rate 
success_1=[ 0.956 0.562 0.078 0 0 ];   %success_rate 
success_2=[1 0.999 0.977 0.745 0.266];   %success_rate 
success_4=[ 1 1 0.993 0.927 0.603  ];   %success_rate 
 %1 0.98 0.87 0.56 0.12 

%figure;
subplot(1,4,2)
h1=plot(B,success_opt,B,success_1,B,success_2,B,success_4);
title('(b) Fat-tree, N=12');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'Optimal','LB(K=1)','LB(K=2)','LB(K=4)','SouthWest');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 


%BCube
N=[2,4,6,8,9,10];%B=100-700;
success_opt=[ 1 1 1 1 1 1 ];   %success_rate 
success_1=[1 1 1 0.66 0.37  0.161 ];   %success_rate 
success_2=[  1 1 1 0.97 0.88  0.696  ];   %success_rate 
success_4=[  1 1 1 0.99 0.93 0.743 ];   %success_rate 
 %1 0.98 0.87 0.56 0.12 

%figure;
subplot(1,4,3);
h1=plot(N,success_opt,N,success_1,N,success_2,N,success_4);
title('(c) BCube, B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'Optimal','LB(K=1)','LB(K=2)','LB(K=4)','SouthWest');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%BCube
N=12;
B=[200,250,300,350,400];
success_opt=[  1 1 1 1 0.999  ];   %success_rate 
success_1=[0.988 0.754 0.114 0.003 0 ];   %success_rate 
success_2=[ 1 1 0.898 0.449 0.04 ];   %success_rate 
success_4=[ 0.999 0.999 0.964 0.533 0.071  ];   %success_rate 
 %1 0.98 0.87 0.56 0.12 

%figure;
subplot(1,4,4);

h1=plot(B,success_opt,B,success_1,B,success_2,B,success_4);
title('(d) BCube, N=12');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'Optimal','LB(K=1)','LB(K=2)','LB(K=4)','SouthWest');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 


