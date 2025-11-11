
clear;clc;close all
coeff=[ 1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];


kh=linspace(0,2*pi,100);

(2*(2^(1/3)+2^(2/3))/(sqrt(2)*2*sum(abs(coeff))))
2/(sqrt(2)*2*sum(abs(coeff)))
sqrt(3)/(sqrt(2)*sum(abs(coeff)))


   % 1.4455
    % 0.8793


coeff = [1.55427, -0.27162, 0.07394, -0.0208847, 0.005117, -0.000919339, 0.0000882575];

(2*(2^(1/3)+2^(2/3))/(sqrt(2)*2*sqrt(sum(abs(coeff)))))


sqrt(3)/(sqrt(2)*sqrt(sum(abs(coeff))))


  % 1.4504
    % 0.8823   slighterly better 
