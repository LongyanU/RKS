clear;
clc;
close all

M=7;

kh=linspace(0.01,0.7*pi,M);

% kh=linspace(0.01,0.6*pi,M);
% M=9;
kh=linspace(0.001,0.6*pi,M);
AA=zeros(M,M);
b=zeros(M,1);

for ii=1:M %kx=1
    
    for kk=1:M
        AA(ii,kk)=2*(sin((kk-0.5)*kh(ii)))*2*(sin(0.5*kh(ii)))*(2*cos(kh(ii))-2);
    end
  
    b(ii)=-(kh(ii))^4;
end

c=AA\b;%ÇóÏµÊý
length(c)
% temp1=-2*sum(c);
digits(6)
vpa(c)'


kh=linspace((pi)/(100),(pi),100);

a=0;
for m=1:length(c)
    a=a+2*c(m)*(sin((m-.5)*kh))*2.*(sin(0.5*kh)).*(2*cos(kh)-2);
end

hold on;plot(-a-kh.^4,'k','LineWidth',2);grid on
axis([0 100 -2.5*10^-4 5.5*10^-4])
legend('The previous SGFD method','The proposed SGFD method')
xlabel('The percentage of kh')
ylabel('Error')
