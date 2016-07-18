clear all;
close all;


%BCube
N=[2,4,6,8,9,10];%B=100-700;
success_opt=[ 1 1 1 1 1 1 ];   %success_rate 
success_1=[1 1 1 0.66 0.37  0.161 ];   %success_rate 
success_2=[  1 1 1 0.97 0.88  0.696  ];   %success_rate 
success_4=[  1 1 1 0.99 0.93 0.743 ];   %success_rate 
 %1 0.98 0.87 0.56 0.12 

figure;
h1=plot(N,success_opt,N,success_1,N,success_2,N,success_4);
title('BCube, B=100-700');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'Optimal','Load-balance(K=1)','Load-balance(K=2)','Load-balance(K=4)');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',1.5);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%BCube
N=[2,4,6,8,9,10]';%B=100-700;
runtime_opt=[  2.595 7.265 13.7957 25.8783 31.994 41.9752 ]';   %runtime_rate 
runtime_4=[ 0.025 0.413636 5.29098 10.1791 11.8209 13.3308 ]';   %runtime_rate 
runtime_2=[ 0.02 0.075 1.7238 5.38728 6.08834 6.95003 ]';   %runtime_rate 
runtime_1=[0.0197044 0.0413534 0.3447 0.700005 0.880557 0.836443 ]';   %runtime_rate 

 %1 0.98 0.87 0.56 0.12 

figure;
h1=plot(N,runtime_opt,N,runtime_1,N,runtime_2,N,runtime_4);
title('BCube, B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Running time of routing algorith(ms)');
legend(h1,'Optimal','Load-balance(K=1)','Load-balance(K=2)','Load-balance(K=4)');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',1.5);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 


%Fat-tree
N=[2,4,6,8,9,10];%B=100-700Mbps;
success_opt=[  1 1 1 1 0.995 0.99  ];   %success_rate 
success_1=[1 1 0.84 0.42 0.23 0.085 ];   %success_rate 
success_2=[  1 1 0.995 0.96 0.925 0.75  ];   %success_rate 
success_4=[ 1 1 1 0.99 0.975 0.905  ];   %success_rate 
 %1 0.98 0.87 0.56 0.12 

figure;
h1=plot(N,success_opt,N,success_1,N,success_2,N,success_4);
title('Fat-tree, B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'Optimal','Load-balance(K=1)','Load-balance(K=2)','Load-balance(K=4)');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',1.5);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

