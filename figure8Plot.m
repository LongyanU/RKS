clear;
close all;
clc

% clip level 0.5
load('figure8a_TraRK.mat')
plotimage((1:nt-2)*2.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recorda=seis_record;



load('figure8b_Non_balance_intime_2.mat')
plotimage((1:nt-2)*2.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordb=seis_record;


load('figure8c_RK_non_balance_inspace.mat')
plotimage((1:nt-2)*2.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordc=seis_record;



load('figure8d_nonbalance_TS.mat')
plotimage((1:nt-2)*2.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordd=seis_record;



figure;plot((1:nt-2)*2.7,seis_recorda((1:nt-2),75),'c','linewidth',2)

hold on;plot((1:nt-2)*2.7,seis_recordb((1:nt-2),75),'r','linewidth',2)
hold on;plot((1:nt-2)*2.7,seis_recordc((1:nt-2),75),'k','linewidth',2)

hold on;plot((1:nt-2)*2.7,seis_recordd((1:nt-2),75),'g','linewidth',2)



xlabel('time(ms)')
ylabel('Amp')
% legend('traditional FD withoutRK','Tra RK','RK_NonBalanceInSpace','RK_NonBalanceInBoth','TDE')
legend('Tra RK','RK-NonBalanceInTime','RK-NonBalanceInSpace','RK-NonBalanceInBoth')
grid on
box on


axis([ 0 5300 -14 21])