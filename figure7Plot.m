clear;
close all;
clc

% clip level 0.5


load('figure7a_TraRK.mat')
plotimage((1:nt-2)*1.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recorda=seis_record;



load('figure7b_Non_balance_intime_2.mat')
plotimage((1:nt-2)*1.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordb=seis_record;


load('figure7c_RK_non_balance_inspace.mat')
plotimage((1:nt-2)*1.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordc=seis_record;



load('figure7d_nonbalance_TS.mat')
plotimage((1:nt-2)*1.5,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordd=seis_record;



figure;plot((1:nt-2)*1.5,seis_recorda((1:nt-2),75),'c','linewidth',2)

hold on;plot((1:nt-2)*1.5,seis_recordb((1:nt-2),75),'r','linewidth',2)
hold on;plot((1:nt-2)*1.5,seis_recordc((1:nt-2),75),'k','linewidth',2)

hold on;plot((1:nt-2)*1.5,seis_recordd((1:nt-2),75),'g','linewidth',2)


load('figure7e_Tra.mat')
hold on;plot((1:nt-2)*1.5,seis_record((1:nt-2),75),'b','linewidth',2)

load('figure7f_TDE.mat')
NT=nt;
B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix

temp=I*seis_record(:,75);

hold on;plot((1:nt-2)*1.5,temp((1:nt-2)),'m','linewidth',2);


xlabel('time(ms)')
ylabel('Amp')
% legend('traditional FD withoutRK','Tra RK','RK_NonBalanceInSpace','RK_NonBalanceInBoth','TDE')
legend('Tra RK','RK-NonBalanceInTime','RK-NonBalanceInSpace','RK-NonBalanceInBoth','without RK','TDE')
grid on
box on
% axis([ 2190 2250 -8*10^-3 6*10^-3])

axis([ 0 5300 -14 21])