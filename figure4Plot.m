clear;
close all;
clc

% clip level 0.5

load('figure4a.mat')

plotimage((1:nt-2)*3.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recorda=seis_record;



load('figure4b.mat')
% load('figure2b_single_precesion.mat')

% load('figure2a_FreeSurface_2.mat')
plotimage((1:nt-2)*3.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordb=seis_record;


% load('figure2b_2.mat')
% % load('figure2a_FreeSurface_2.mat')
% plotimage((1:nt-2)*3.7,seis_record(1:nt-2,50:end-50))
% xlabel('x/dx')
% ylabel('Time(ms)')
% title('')
% seis_recorde=seis_record;


load('figure4c.mat')
plotimage((1:nt-2)*3.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordc=seis_record;


% load('figure2d_FreeSurface.mat')

load('figure4d.mat')
plotimage((1:nt-2)*3.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordd=seis_record;



figure;plot((1:nt-2)*3.7,seis_recorda((1:nt-2),75),'c','linewidth',2)

hold on;plot((1:nt-2)*3.7,seis_recordb((1:nt-2),75),'r','linewidth',2)
hold on;plot((1:nt-2)*3.7,seis_recordc((1:nt-2),75),'k','linewidth',2)

hold on;plot((1:nt-2)*3.7,seis_recordd((1:nt-2),75),'g--','linewidth',2)

% hold on;plot((1:nt-2)*3.7,seis_recorde((1:nt-2),75),'y--','linewidth',2)



xlabel('time(ms)')
ylabel('Amp')
% legend('traditional FD withoutRK','Tra RK','RK_NonBalanceInSpace','RK_NonBalanceInBoth','TDE')
% legend('Tra RK','RK-NonBalanceInTime','RK-NonBalanceInSpace','RK-NonBalanceInBoth','without RK','TDE')
% legend('Tra RK','RK-NonBalanceInTime','RK-NonBalanceInSpace','RK-NonBalanceInBoth','without RK','TDE','LW')
legend('Tra RKS','NU-T-RKS','NU-S-RKS','NU-TS-RKS')
grid on
box on
% axis([ 2190 2250 -8*10^-3 6*10^-3])

axis([ 0 1995*2.1 -30 30])