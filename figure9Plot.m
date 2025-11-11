clear;
close all;
clc


% clip level 0.5
 load('figure9a_TraRK.mat')
plotimage((1:nt-2)*4.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recorda=seis_record;





load('figure9b_RK_non_balance_inspace.mat')
plotimage((1:nt-2)*4.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordc=seis_record;






figure;plot((1:nt-2)*4.5,seis_recorda((1:nt-2),75),'c','linewidth',2)

hold on;plot((1:nt-2)*4.5,seis_recordc((1:nt-2),75),'k','linewidth',2)



xlabel('time(ms)')
ylabel('Amp')
% legend('traditional FD withoutRK','Tra RK','RK_NonBalanceInSpace','RK_NonBalanceInBoth','TDE')
legend('Tra RK','RK-NonBalanceInSpace')
grid on
box on
% axis([ 2190 2250 -8*10^-3 6*10^-3])

axis([ 0 5300 -11.5 19])