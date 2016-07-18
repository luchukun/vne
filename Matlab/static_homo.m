clear all;
close all;


%FatTree
N=[2,4,6,8,9,10];
suc_hvc=[1 1 0.68 0.57 0.435 0.43  ];   %suc_rate HVC_ACE
suc_pert_opt=[1 1 0.96 0.89 0.825 0.76 ];   %suc_rate fisrt fit
figure;
subplot(2,2,1)
h1=plot(N,suc_hvc,N,suc_pert_opt);
title('(a) Fat-tree, B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'HVC-ACE','pertubation','Location','SouthWest');
set(gca,'XLim',[2,  10]);set(gca,'XGrid','on','YGrid','on');
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'},{'MarkerSize'},{10;8});
set(h1,'LineWidth',2); 
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%FatTree
B=[200 250 300 350 400] ;
suc_hvc=[0.995 0.985 0.955 0.825 0.735  ];   %suc_rate HVC_ACE
suc_pert_opt=[1 1 1 1 0.99 ];   %suc_rate fisrt fit
%figure;
subplot(2,2,3)
h1=plot(B,suc_hvc,B,suc_pert_opt);
title('(c) Fat-tree, N=8');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'HVC-ACE','pertubation','Location','SouthWest');
set(gca,'XLim',[200,400]);set(gca,'XGrid','on','YGrid','on');
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'},{'MarkerSize'},{10;8});
set(h1,'LineWidth',2); 
%Bcube
N=[2,4,6,8,9,10];
suc_hvc=[1 1 0.915 0.68 0.62 0.59];   %suc_rate HVC_ACE
suc_pert_opt=[1 1 1 0.96 0.91 0.855 ];   %suc_rate fisrt fit
%figure;
subplot(2,2,2)
h1=plot(N,suc_hvc,N,suc_pert_opt);
title('(b) BCube, B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'HVC-ACE','pertubation','Location','SouthWest');
set(gca,'XLim',[2,  10]);set(gca,'XGrid','on','YGrid','on');
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'},{'MarkerSize'},{10;8});
set(h1,'LineWidth',2); 
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%Bcube
B=[200 250 300 350 400] ;
suc_hvc=[1 1 0.985 0.945 0.885  ];   %suc_rate HVC_ACE
suc_pert_opt=[ 1 1 1 1 1 ];   %suc_rate fisrt fit
%figure;
subplot(2,2,4)
h1=plot(B,suc_hvc,B,suc_pert_opt);
title('(d) BCube, N=8');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'HVC-ACE','pertubation','Location','SouthWest');
set(gca,'XLim',[200,400]);set(gca,'XGrid','on','YGrid','on');
set(h1,{'Color'},{'r';'b'});
set(h1,{'Marker'},{'*';'O'},{'MarkerSize'},{10;8});
set(h1,'LineWidth',2); 
