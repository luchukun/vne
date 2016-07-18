clear all;
close all;



%Bcube
B= [200 250 300 350 400 ]

thread1=[  21.679 34.2959 44.8063 53.8786 64.4589  ];   
thread2=[16.1176 25.3533 32.784 40.7359 49.6098 ];   
thread3=[13.8469 22.1727 27.303 35.279 38.7603 ];   
thread4=[ 12.6016 19.4198 23.9291 27.7123 32.9672 ];   
thread5=[12.9635 18.6558 25.7094 30.3546 32.2756 ];   
thread6=[13.5993 19.7836 23.9182 29.0601 32.5649 ];   
thread7=[ 14.1275 20.1899 23.959 28.4524 31.2683  ];   
thread8=[13.7336 21.7929 25.3233 31.4487 34.4515 ];   

figure;
h1=plot(B,thread1,B,thread2,B,thread3,B,thread4,B,thread5,B,thread6,B,thread7,B,thread8);
title('BCube,N=12');

xlabel('Average bandwidth demand');
ylabel('Running time(ms)');
legend(h1,'thread=1','thread=2','thread=3','thread=4','thread=5','thread=6','thread=7','thread=8','Location','SouthWest');
set(gca,'XLim',[200,400]);set(gca,'XGrid','on','YGrid','on')
%set(h1,{'Color'},{'m';'g';'c';'b'});
%set(h1,{'Marker'},{'s';'+';'x';'O'},{'MarkerSize'},{8;10;10;8});
set(h1,'LineWidth',2);
%set(h1,'LineWidth',1.5,{'LineStyle'},{'-';'-';'-';}) 

