clear all;
close all;


%Fat-tree
B=[200,250,300,350,400];
suc_randrop=[ 0.999 0.989 0.949 0.711 0.216 ];   %success_rate backtracking
suc_pert=[ 1 0.999 0.977 0.745 0.266 ];   %success_rate pertuabtion 

figure;
h1=plot(B,suc_randrop,B,suc_pert);
title('Fat-tree');

xlabel('Average bandwidth deband');
ylabel('Success rate');
legend(h1,'ramdom drop','pertubation');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'});
set(h1,'LineWidth',1);

%BCube
B=[200,250,300,350,400];
suc_randrop=[ 0.999 0.994 0.858 0.353 0.038 ];   %success_rate backtracking
suc_pert=[1 1 0.898 0.449 0.04 ];   %success_rate pertuabtion 

figure;
h1=plot(B,suc_randrop,B,suc_pert);
title('BCube,K=2');

xlabel('Average bandwidth deband');
ylabel('Success rate');
legend(h1,'ramdom drop','pertubation');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'});
set(h1,'LineWidth',1);

%Vl2
B=[200,250,300,350,400];
suc_randrop=[0.997 0.991 0.96 0.765 0.338  ];   %success_rate backtracking
suc_pert=[1 0.997 0.974 0.834 0.384 ];   %success_rate pertuabtion 

figure;
h1=plot(B,suc_randrop,B,suc_pert);
title('VL2,K=4');

xlabel('Average bandwidth deband');
ylabel('Success rate');
legend(h1,'ramdom drop','pertubation');
set(gca,'XLim',[200,400]);
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'});
set(h1,'LineWidth',1);


