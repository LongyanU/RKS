clear;
clc;
close all

d= [1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];

kh=linspace(1/1000,pi,100);

cdcfinal=zeros(5,100);
D=0;


h=10;
v=2500;
% v=5000;
v=4900;
tau=0.001
dt=tau;
r=v*tau/h;

for ii=1:100

    for aa=1:5

        xita=(aa-1)*pi/16;
        tempyy=100;
        for cdc=0.79:0.00001:1.2
            tempD1=0;
            for n=1:7
                tempD1=tempD1+2*d(n)* sin((n-0.5)*kh(ii)*sin(xita));
            end

            tempD2=0;
            for n=1:7
                tempD2=tempD2+2*d(n)*sin((n-0.5)*kh(ii)*cos(xita));
            end

            D=sqrt(tempD1^2+tempD2^2);

            temp1=2*cos(kh(ii)*r*cdc) -2+ (r*D  - (r*D)^3/12)*(r*D);

            if abs(temp1)<tempyy
                tempyy=abs(temp1);
                cdcfinal(aa,ii)=cdc;
            end

        end
    end
end


% a1=real(h/v*(1./(cdcfinal)-1));

a1=cdcfinal;
for ii=1:5
    if ii==1
        figure;plot(100*kh/(pi),a1(ii,:),'r','LineWidth',1)
    elseif ii==2
        hold on;plot(100*kh/(pi),a1(ii,:),'b','LineWidth',1)
    elseif ii==3
        hold on;plot(100*kh/(pi),a1(ii,:),'c','LineWidth',1)
    elseif ii==4
        hold on;plot(100*kh/(pi),a1(ii,:),'k','LineWidth',1)
    else
        hold on;plot(100*kh/(pi),a1(ii,:),'m','LineWidth',1)
    end
end

grid on
legend('\theta=0','\theta=\pi/16','\theta=2\pi/16','\theta=3\pi/16','\theta=4\pi/16')

xlabel('percentage of kh')
ylabel('$\frac{v_{\mathrm{fd}}}{v}$', 'Interpreter', 'latex');
% xlim([0,pi])
% ylim([0.89,1.02])
% ylim([0.94,1.08])
% axis([0 pi -10*10^-5 3*10^-5])

% axis([0 100 -5*10^-5 5*10^-5])
% axis([0 100 -4*10^-5 7*10^-5])
% axis([0 100 -0.5*10^-5 1*10^-4])

axis([0 100 0.9 1.14])