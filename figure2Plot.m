clear;
close all;
clc

% clip level 0.5









 load('figure2a.mat')
plotimage((1:nt-2)*2.1,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recorda=seis_record;



load('figure2b.mat')
plotimage((1:nt-2)*2.1,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordb=seis_record;


load('figure2c.mat')
plotimage((1:nt-2)*2.1,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordc=seis_record;


% load('figure2d.mat')
load('figure2d_2.mat')
plotimage((1:nt-2)*2.1,seis_record(1:nt-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordd=seis_record;


figure;plot((1:nt-2)*2.1,seis_recorda((1:nt-2),75),'c','linewidth',2)

hold on;plot((1:nt-2)*2.1,seis_recordb((1:nt-2),75),'r','linewidth',2)
hold on;plot((1:nt-2)*2.1,seis_recordc((1:nt-2),75),'k','linewidth',2)

hold on;plot((1:nt-2)*2.1,seis_recordd((1:nt-2),75),'g','linewidth',2)


load('figure2e_2order_SGFD.mat')
hold on;plot((1:nt-2)*2.1,seis_record((1:nt-2),75),'b','linewidth',2)

load('figure2f_TDE.mat')
NT=nt;
B = ifft(exp(-2i*sin([0:NT-1]*pi/(2*NT))'*[0:NT-1]).*cos([0:NT-1]'*pi/(2*NT)),2*NT,'symmetric');
I = B(1:NT,1:NT)'; % <- The Inverse Time Dispersion Transform matrix

temp=I*seis_record(:,75);

hold on;plot((1:nt-2)*2.1,temp((1:nt-2)),'m','linewidth',2);


xlabel('time(ms)')
ylabel('Amp')
% legend('traditional FD withoutRK','Tra RK','RK_NonBalanceInSpace','RK_NonBalanceInBoth','TDE')
legend('Tra RK','RK-NonBalanceInTime','RK-NonBalanceInSpace','RK-NonBalanceInBoth','without RK','TDE')
grid on
box on
% axis([ 2190 2250 -8*10^-3 6*10^-3])

axis([ 0 3130 -18 24.5])