
clear;clc
close all

% Tra R_K [22/24 1/24 1/24]
A=[1 1 1;
    3 1 -1;
    9  1 1];

rank(A)
cond(A)

b=[1 1 4/3];

c=A\b'

digits(6)
vpa(c')
vpa([22/24 1/24 1/24])