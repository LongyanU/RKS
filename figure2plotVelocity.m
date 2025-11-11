clear;clc;close all

% load('figure3_tra.mat')

load('figure2a.mat')

figure;imagesc(v)
xlabel('x/dx')
ylabel('z/dz')

h = colorbar;
set(get(h,'title'),'string','m/s');
hold on ;plot(nx/2,15+45,'*r')
% annotation('textbox', [0.93 0.89 0.05 0.05], ...
%            'String', 'm/s');