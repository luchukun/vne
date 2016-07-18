clear all;
close all;


%VL2
N=12;
B=[200,250,300,350,400];
mlu_opt=[0.955069 0.976387 0.979974 0.977653 0.981051];   %mlu_rate 
mlu_1=[0.975146 0.983558 0.986838 0.988516 0.992392 ];   %mlu_rate 
mlu_2=[0.966939 0.982205 0.983874 0.984812 0.987772 ];   %mlu_rate 
mlu_4=[0.970678 0.983016 0.983584 0.983853 0.987613 ];   %mlu_rate 
mlu_8=[0.962075 0.977817 0.980775 0.980947 0.983869 ];   %mlu_rate
 %1 0.98 0.87 0.56 0.12 

figure;
h1=plot(B,mlu_opt,B,mlu_1,B,mlu_2,B,mlu_4,B,mlu_8);
title('VL2, N=12');

xlabel('Average bandwidth deband');
ylabel('MLU');
legend(h1,'Optimal','Load-balance(K=1)','Load-balance(K=2)','Load-balance(K=4)','Load-balance(K=8)');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'k';'m'});
set(h1,{'Marker'},{'*';'O';'s';'+';'x'});
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%Bcube
N=12;
B=[200,250,300,350,400];
mlu_opt=[ 0.64171 0.779884 0.878025 0.943125 0.975616 ];   %mlu_rate 
mlu_1=[0.988 0.754 0.114 0.003 0 ];   %mlu_rate 
mlu_2=[0.898059 0.952899 0.973275 0.98121 0.980187 ];   %mlu_rate 
mlu_4=[ 0.931909 0.976312 0.981183 0.984365 0.981416 ];   %mlu_rate 
 %1 0.98 0.87 0.56 0.12 

figure;
h1=plot(B,mlu_opt,B,mlu_1,B,mlu_2,B,mlu_4);
title('Bcube, N=12');

xlabel('Average bandwidth deband');
ylabel('MLU');
legend(h1,'Optimal','Load-balance(K=1)','Load-balance(K=2)','Load-balance(K=4)');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'});
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 


%FatTree
N=12;
B=[200,250,300,350,400];
mlu_opt=[1 0.998 0.999 0.993 0.983 ];   %mlu_rate 
mlu_1=[0.979685 0.986444 0.991088 0 0 ];   %mlu_rate 
mlu_2=[ 0.967794 0.979343 0.9833 0.985825 0.988804 ];   %mlu_rate 
mlu_4=[ 0.960526 0.979796 0.983901 0.985245 0.984827  ];   %mlu_rate 
 %1 0.98 0.87 0.56 0.12 

figure;
h1=plot(B,mlu_opt,B,mlu_1,B,mlu_2,B,mlu_4);
title('FatTree, N=12');

xlabel('Average bandwidth deband');
ylabel('MLU');
legend(h1,'Optimal','Load-balance(K=1)','Load-balance(K=2)','Load-balance(K=4)');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b';'g';'k'});
set(h1,{'Marker'},{'*';'O';'s';'+'});
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 
