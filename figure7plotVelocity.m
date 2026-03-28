clear;clc;close all

% load('figure3_tra.mat')

load('figure1_TraRK.mat')

figure;imagesc(c1)
xlabel('x/dx')
ylabel('z/dz')

h = colorbar;
set(get(h,'title'),'string','m/s');
hold on ;plot(nx/2,1,'*r')
% annotation('textbox', [0.93 0.89 0.05 0.05], ...
%            'String', 'm/s');