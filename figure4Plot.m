clear;
close all;
clc

% clip level 0.5

load('figure4a_TraRK.mat')
plotimage((1:nt-2)*3.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recorda=seis_record;



load('figure4b_Non_balance_intime.mat')
% load('figure1_Non_balance_intime.mat')
plotimage((1:nt-2)*3.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordb=seis_record;


load('figure4c_RK_non_balance_inspace.mat')
plotimage((1:nt-2)*3.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordc=seis_record;



load('figure4d_nonbalance_TS.mat')
plotimage((1:nt-2)*3.7,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordd=seis_record;



figure;plot((1:nt-2)*3.7,seis_recorda((1:nt-2),75),'c','linewidth',2)

hold on;plot((1:nt-2)*3.7,seis_recordb((1:nt-2),75),'r','linewidth',2)
hold on;plot((1:nt-2)*3.7,seis_recordc((1:nt-2),75),'k','linewidth',2)

hold on;plot((1:nt-2)*3.7,seis_recordd((1:nt-2),75),'g','linewidth',2)


xlabel('time(ms)')
ylabel('Amp')
legend('Tra RK','RK-NonBalanceInTime','RK-NonBalanceInSpace','RK-NonBalanceInBoth')
grid on
box on

axis([ 0 3130 -13.8 26.9])