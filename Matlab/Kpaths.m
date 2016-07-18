clear all;
close all;


%VL2
N=12;
B=[200,250,300,350,400];
success_opt=[ 1 1 1 1 0.97 ];   %success_rate 
success_1=[1 0.982 0.781 0.313 0.035 ];   %success_rate 
success_2=[0.998 0.995 0.946 0.682 0.275 ];   %success_rate 
success_4=[1 0.997 0.974 0.834 0.384 ];   %success_rate 
success_8=[1 0.999 0.986 0.918 0.602 ];   %success_rate
 %1 0.98 0.87 0.56 0.12 

figure;
h1=plot(B,success_opt,B,success_1,B,success_2,B,success_4,B,success_8);
title('VL2, N=12');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'Optimal','Load-balance(K=1)','Load-balance(K=2)','Load-balance(K=4)','Load-balance(K=8)');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'k';'m'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'},{'MarkerSize'},{10;8;8;10;10});
set(h1,'LineWidth',1.5);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%BCube
N=12;
B=[200,250,300,350,400];
success_opt=[  1 1 1 1 0.999  ];   %success_rate 
success_1=[0.988 0.754 0.114 0.003 0 ];   %success_rate 
success_2=[ 1 1 0.898 0.449 0.04 ];   %success_rate 
success_4=[ 0.999 0.999 0.964 0.533 0.071  ];   %success_rate 
 %1 0.98 0.87 0.56 0.12 

figure;
h1=plot(B,success_opt,B,success_1,B,success_2,B,success_4);
title('BCube, N=12');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'Optimal','Load-balance(K=1)','Load-balance(K=2)','Load-balance(K=4)');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',1.5);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 


%Fat-tree
N=12;
B=[200,250,300,350,400];
success_opt=[1 0.998 0.999 0.993 0.983 ];   %success_rate 
success_1=[ 0.956 0.562 0.078 0 0 ];   %success_rate 
success_2=[1 0.999 0.977 0.745 0.266];   %success_rate 
success_4=[ 1 1 0.993 0.927 0.603  ];   %success_rate 
 %1 0.98 0.87 0.56 0.12 

figure;
h1=plot(B,success_opt,B,success_1,B,success_2,B,success_4);
title('Fat-tree, N=12');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'Optimal','Load-balance(K=1)','Load-balance(K=2)','Load-balance(K=4)');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'},{'MarkerSize'},{10;8;8;10;});
set(h1,'LineWidth',1.5);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 
