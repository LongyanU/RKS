clear;
close all;
clc

% clip level 0.5




load('figure5a_TraRK.mat')
plotimage((1:nt-2)*6.1,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recorda=seis_record;


load('figure5b_non_balance_inspace_2.mat')
plotimage((1:nt-2)*6.1,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordc=seis_record;



figure;plot((1:nt-2)*6.1,seis_recorda((1:nt-2),75),'c','linewidth',2)

hold on;plot((1:nt-2)*6.1,seis_recordc((1:nt-2),75),'k','linewidth',2)




xlabel('time(ms)')
ylabel('Amp')

legend('Tra RK','RK-NonBalanceInSpace')
grid on
box on


axis([ 0 3100 -13.8 28.5])