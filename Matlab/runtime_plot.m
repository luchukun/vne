
N=[2,4,6,8,10];
success_b=[  264 516 732 11556 188650 517344 ];   %success_rate backtracking
success_p=[  210 430 610 840 1650 2900 ];   %success_rate pertuabtion
success_f=[  220 420 550 580 550 540];   %success_rate fisrt fit
 
figure;
h1=plot(N,success_b,N,success_p,N,success_f);
title('BCube, K=2,B=100~700');

xlabel('Number of VMs (N)');
ylabel('Running time/ms');
legend(h1,'backtracking','pertubation','first-fit');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';});
set(h1,{'Marker'},{'*';'O';'s';});

%running time
N=[2,4,6,8,10];
success_b=[ 68 130 402 3662 54818];   %success_rate backtracking
success_p=[ 64 128 190 320 586 ];   %success_rate pertuabtion
success_f=[78 158 162 174 150];   %success_rate fisrt fit
figure;
h1=plot(N,success_b,N,success_p,N,success_f);
title('Fat-tree, K=2,B=100~700');

xlabel('Number of VMs (N)');
ylabel('Running time/ms');
legend(h1,'backtracking','pertubation','first-fit');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';});
set(h1,{'Marker'},{'*';'O';'s';});



N=[2,4,6,8,10];
success_b=[ 52 126 206 1572 27274 ];   %success_rate backtracking
success_p=[ 54 106 172 264 420  ];   %success_rate pertuabtion
success_f=[ 64 134 136 122 138 ];   %success_rate fisrt fit
figure;
h1=plot(N,success_b,N,success_p,N,success_f);
title(VL2, K=2,B=100~700');

xlabel('Number of VMs (N)');
ylabel('Running time/');
legend(h1,'backtracking','pertubation','first-fit');
set(gca,'XLim',[2,10]);
set(h1,{'Color'},{'r';'b';'g';});
set(h1,{'Marker'},{'*';'O';'s';});