clear;
close all;
clc

% clip level 0.5

load('figure2a_FreeSurface_2.mat')

plotimage((1:2000-2)*2.1,seis_record(1:2000-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recorda=seis_record;



% load('figure2b_FreeSurface.mat')
load('figure2b_2.mat')
% load('figure2a_FreeSurface_2.mat')
plotimage((1:2000-2)*2.1,seis_record(1:2000-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')
seis_recordb=seis_record;


load('figure2c.mat')
plotimage((1:2000-2)*2.1,seis_record(1:2000-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordc=seis_record;


% load('figure2d_FreeSurface.mat')

load('figure2d.mat')
% load('figure2d_2.mat')
plotimage((1:nt-2)*2.1,seis_record(1:2000-2,50:end-50))
xlabel('x/dx')
ylabel('Time(ms)')
title('')

seis_recordd=seis_record;


figure;plot((1:nt-2)*2.1,seis_recorda((1:nt-2),75),'c','linewidth',2)

hold on;plot((1:nt-2)*2.1,seis_recordb((1:nt-2),75),'r-.','linewidth',2)
hold on;plot((1:nt-2)*2.1,seis_recordc((1:nt-2),75),'k','linewidth',2)

hold on;plot((1:nt-2)*2.1,seis_recordd((1:nt-2),75),'g','linewidth',2)


load('figure2e_2order_SGFD.mat')
hold on;plot((1:nt-2)*2.1,seis_record((1:nt-2),75),'b','linewidth',2)

%% --- Corrected TDE Inverse Transform Block ---
load('figure2f.mat')
% 1. Get the actual vector and its length
trace = seis_record(:, 75);
L = length(trace);      % Total number of time samples
K = L - 1;              % Maximum index for the [0:K] range

% 2. Reconstruct the matrix I to be exactly L x L
% Use (0:K) to ensure we have exactly L elements
idx = (0:K)'; 
% Matrix construction (ensure dimensions are L x L)
B = ifft(exp(-2i * sin(idx * pi / (2 * K)) * idx') .* cos(idx * pi / (2 * K)), 2 * K, 'symmetric');
I_mat = B(1:L, 1:L)'; 

% 3. Perform the multiplication
% Ensure trace is a column vector
if size(trace, 2) > size(trace, 1)
    trace = trace'; 
end

temp = I_mat * trace;

% 4. Plot
hold on;
plot((1:nt-2)*2.1, temp(1:nt-2), 'm--', 'linewidth', 2);



load('figure2fabs.mat')
hold on;plot((1:nt-2)*1.4,seis_record((1:nt-2),75),'y--','linewidth',2)

% load('figure2fabs_2.mat')
% hold on;plot((1:nt-2)*1.4,seis_record((1:nt-2),75),'c--','linewidth',2)

%  load('figure2d_rbig.mat')
% hold on;plot((1:nt-2)*3.7,seis_record((1:nt-2),75),'b','linewidth',2)

xlabel('time(ms)')
ylabel('Amp')
% legend('traditional FD withoutRK','Tra RK','RK_NonBalanceInSpace','RK_NonBalanceInBoth','TDE')
% legend('Tra RK','RK-NonBalanceInTime','RK-NonBalanceInSpace','RK-NonBalanceInBoth','without RK','TDE')
% legend('Tra RK','RK-NonBalanceInTime','RK-NonBalanceInSpace','RK-NonBalanceInBoth','without RK','TDE','LW')
legend('Tra RKS','NU-T-RKS','NU-S-RKS','NU-TS-RKS','without RKS','TDE','Tra ABS')
grid on
box on
% axis([ 2190 2250 -8*10^-3 6*10^-3])

axis([ 0 1995*2.1 -30 30])