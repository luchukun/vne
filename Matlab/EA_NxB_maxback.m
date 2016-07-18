clear all;
close all;



%Bcube
N=[2,4,6,8,9,10];
%suc_bktrack=[  1 1 1 0.99 0.96 0.82   ];   %success_rate backtracking
suc_pert=[1 1 0.999 0.986 0.92 0.744 ];   %success_rate pertuabtion
suc_pert1n=[  1 1 0.991 0.889 0.763 0.528  ];   %success_rate fisrt fit
suc_pert2n=[1 1 0.999 0.941 0.848 0.634  ];   %success_rate fisrt fit
suc_pert3n=[1 1 0.999 0.982 0.877 0.674 ];   %success_rate fisrt fit 
figure;
subplot(1,2,1);
h1=plot(N,suc_pert1n,N,suc_pert2n,N,suc_pert3n,N,suc_pert);
title('(a) BCube,B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Success rate');
legend(h1,'\Delta=N','\Delta=2N','\Delta=3N','\Delta=(|V_s|-1)N','Location','SouthWest');
set(gca,'XLim',[2,10]);set(gca,'XGrid','on','YGrid','on')
set(h1,{'Color'},{'m';'g';'c';'b'});
set(h1,{'Marker'},{'s';'+';'x';'O'},{'MarkerSize'},{8;10;10;8});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%

%Bcube
N=8;
B=[200 250 300 350 400] ;
%suc_bktrack=[  1 1 1 1 0.99 ];   %success_rate backtracking
suc_pert=[1 1 1 0.997 0.986  ];   %success_rate pertuabtion
suc_pert1n=[   0.998 0.991 0.985 0.934 0.808 ];   %success_rate fisrt fit
suc_pert2n=[1 0.997 0.991 0.968 0.907  ];   %success_rate fisrt fit
suc_pert3n=[  1 1 0.998 0.982 0.946   ];   %success_rate fisrt fit 
%figure;
subplot(1,2,2);
h1=plot(B,suc_pert1n,B,suc_pert2n,B,suc_pert3n,B,suc_pert);
title('(b) BCube,N=8');

xlabel('Average bandwidth demand');
ylabel('Success rate');
legend(h1,'\Delta=N','\Delta=2N','\Delta=3N','\Delta=(|V_s|-1)N','Location','SouthWest');
set(gca,'XLim',[200,400]);set(gca,'XGrid','on','YGrid','on')
set(h1,{'Color'},{'m';'g';'c';'b'});
set(h1,{'Marker'},{'s';'+';'x';'O'},{'MarkerSize'},{8;10;10;8});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 


%Bcube
N=[2,4,6,8,9,10];
%time_bktrack=[  1 1 1 0.99 0.96 0.82   ];   %success_rate backtracking
time_pert=[0.195 0.225 0.56 7.742 14.517 37.404 ];   %success_rate pertuabtion
time_pert1n=[0.194 0.224 0.609 3.943 7.864 15.846  ];   %success_rate fisrt fit
time_pert2n=[0.194 0.224 0.579 5.269 11.597 27.2 ];   %success_rate fisrt fit
time_pert3n=[0.193 0.225 0.557 5.277 15.383 34.901 ];   %success_rate fisrt fit 
figure;
subplot(1,2,1);
h1=plot(N,time_pert1n,N,time_pert2n,N,time_pert3n,N,time_pert);
title('(c) BCube,B=100-700Mbps');

xlabel('Number of VMs (N)');
ylabel('Running time(ms)');
legend(h1,'\Delta=N','\Delta=2N','\Delta=3N','\Delta=(|V_s|-1)N','Location','SouthWest');
set(gca,'XLim',[2,10]);set(gca,'XGrid','on','YGrid','on')
set(h1,{'Color'},{'m';'g';'c';'b'});
set(h1,{'Marker'},{'s';'+';'x';'O'},{'MarkerSize'},{8;10;10;8});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

%

%Bcube
N=8;
B=[200 250 300 350 400] ;
%time_bktrack=[  1 1 1 1 0.99 ];   %success_rate backtracking
time_pert=[ 0.72 1.008 1.719 4.376 11.119 ];   %success_rate pertuabtion
time_pert1n=[ 0.722 1.068 1.51 2.99 5.756  ];   %success_rate fisrt fit
time_pert2n=[ 0.693 0.944 1.631 3.74 7.334] ;   %success_rate fisrt fit
time_pert3n=[ 0.801 0.984 1.613 3.952 8.322  ];   %success_rate fisrt fit 
%figure;
subplot(1,2,2);
h1=plot(B,time_pert1n,B,time_pert2n,B,time_pert3n,B,time_pert);
title('(d) BCube,N=8');

xlabel('Average bandwidth demand');
ylabel('Running time(ms)');
legend(h1,'\Delta=N','\Delta=2N','\Delta=3N','\Delta=(|V_s|-1)N','Location','SouthWest');
set(gca,'XLim',[200,400]);set(gca,'XGrid','on','YGrid','on')
set(h1,{'Color'},{'m';'g';'c';'b'});
set(h1,{'Marker'},{'s';'+';'x';'O'},{'MarkerSize'},{8;10;10;8});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 
