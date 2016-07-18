clear all;
close all;

%Bcube
N=[2,4,6,8,9,10];
suc_hvc=[1 0.996 0.873 0.619 0.541 0.481 ];   %suc_rate HVC_ACE
suc_pert_opt=[1 1 0.991 0.878 0.76 0.696 ];   %suc_rate fisrt fit
figure;
h1=plot(N,suc_hvc,N,suc_pert_opt);
title('BCube');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'HVC-ACE','congestion-aware','Location','SouthWest');
set(gca,'XLim',[2,  10]);
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'});
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%FatTree
N=[2,4,6,8,9,10];
suc_hvc=[1 1 0.689 0.562 0.431 0.41 ];   %suc_rate HVC_ACE
suc_pert_opt=[1 0.99 0.904 0.764 0.668 0.609 ];   %suc_rate fisrt fit
figure;
h1=plot(N,suc_hvc,N,suc_pert_opt);
title('Fat-tree');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'HVC-ACE','congestion-aware','Location','SouthWest');
set(gca,'XLim',[2,  10]);
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'});
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 
