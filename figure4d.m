% 历时 24.829985 秒。
% 历时 26.242479 秒。
clear
clc %%%%%%%
close all
nt=852;
eps=.6;     % stability
isnap=60;    % snapshot sampling

nx=450;
nz=450;

v=ones(nz,nx)*3500;
v(1:nz/2,:)=2200;


p=zeros([nz nx]); pboundarynew=p;pdan=p;
dx=15;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
% dt=0.0038; % calculate time step from stability criterion
dt=0.0037;
tau=dt;
r2=v.*v.*dt*dt/dx/dx;

f0=20*pi;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian


seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=p;
d2pz=p;

% Source location
zs=51;
xs=nz/2;

h=dx;
r=v*dt/h;

p=zeros([nz nx]); Vx=p; Vz=p;
coeff = [1.55427, -0.27162, 0.07394, -0.0208847, 0.005117, -0.000919339, 0.0000882575];

d =[ 0.833333,0.0833333, 0.0833333];
d_prime=d;

dpx_dx=zeros(nz,nx);
dpy_dy=zeros(nz,nx);
Vx2=zeros(nz,nx);
Vz2=zeros(nz,nx);

tic
for it=1:nt-2,

    [dvx_dx, dvy_dy] = Vx_Vz_spatial_derivatives(Vx,Vz,coeff);
    p = p + dt * v.^2.*((dvx_dx + dvy_dy) );

    p(zs,xs)= p(zs,xs)+src(it);
    [p,p]=spongeABC(p,p,nx,nz,50,50,0.009);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p2_prime = p - dt *v.^2.* (dvx_dx + dvy_dy);

    % Stage 3'
    [dpx_dx, dpy_dy] = spatial_derivatives(p);
    vx2_prime = Vx + dt * dpx_dx;
    vy2_prime = Vz + dt * dpy_dy;

    % Stage 4'
    [dvx_dx, dvy_dy] = Vx_Vz_spatial_derivatives(vx2_prime, vy2_prime,coeff);
    % [dvx_dx, dvy_dy] = Vx_Vz_spatial_derivatives_short(vx2_prime, vy2_prime);
    p3_prime = p + dt *v.^2.* (dvx_dx + dvy_dy);

    % Stage 5'
    [dpx_dx1, dpy_dy1] = spatial_derivatives(p2_prime);
    [dpx_dx2, dpy_dy2] = spatial_derivatives(p3_prime);

    % Combine for v_{n+3/2}
    Vx = Vx + dt * (d_prime(1)*dpx_dx + ...
        d_prime(2)*dpx_dx1 + ...
        d_prime(3)*dpx_dx2);

    Vz = Vz + dt * (d_prime(1)*dpy_dy + ...
        d_prime(2)*dpy_dy1 + ...
        d_prime(3)*dpy_dy2);

    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,50,50,0.009);


    if rem(it,60)== 0,
        imagesc(x,z,p), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end
    seis_record(it,:)=p(zs,:);

end
toc
save('figure4d_nonbalance_TS.mat')
function [d1px,d1pz] =spatial_derivatives(p)
h=15;
d1px=(circshift(p,[0 -1])-circshift(p,[0 0]));
d1pz=(circshift(p,[-1])-circshift(p,[0]));
d1px=d1px./h;
d1pz=d1pz./h;
end


function  [d1px,d1pz] =Vx_Vz_spatial_derivatives(Vx,Vz,coeff)
h=15;
d1px11=Vx-circshift(Vx,[0 1]);
d1px12=(circshift(Vx,[0 -1])-circshift(Vx,[0 2]));
d1px13=(circshift(Vx,[0 -2])-circshift(Vx,[0 3]));
d1px14=(circshift(Vx,[0 -3])-circshift(Vx,[0 4]));
d1px15=(circshift(Vx,[0 -4])-circshift(Vx,[0 5]));
d1px16=(circshift(Vx,[0 -5])-circshift(Vx,[0 6]));
d1px17=(circshift(Vx,[0 -6])-circshift(Vx,[0 7]));


d1pz11=Vz-circshift(Vz,[1 0]);
d1pz12=(circshift(Vz,[-1 0])-circshift(Vz,[2 0]));
d1pz13=(circshift(Vz,[-2 0])-circshift(Vz,[3 0]));
d1pz14=(circshift(Vz,[-3 0])-circshift(Vz,[4 0]));
d1pz15=(circshift(Vz,[-4 0])-circshift(Vz,[5 0]));
d1pz16=(circshift(Vz,[-5 0])-circshift(Vz,[6 0]));
d1pz17=(circshift(Vz,[-6 0])-circshift(Vz,[7 0]));


d1px=coeff(1)*d1px11+coeff(2)*d1px12+coeff(3)*d1px13+coeff(4)*d1px14+coeff(5)*d1px15+coeff(6)*d1px16...
    +coeff(7)*d1px17;
d1pz=coeff(1)*d1pz11+coeff(2)*d1pz12+coeff(3)*d1pz13+coeff(4)*d1pz14+coeff(5)*d1pz15+coeff(6)*d1pz16...
    +coeff(7)*d1pz17;
d1px=d1px./h;
d1pz=d1pz./h;
end


function  [d1px,d1pz] =Vx_Vz_spatial_derivatives_short(Vx,Vz)
h=15;
d1px=Vx-circshift(Vx,[0 1]);
d1pz=Vz-circshift(Vz,[1 0]);
d1px=d1px./h;
d1pz=d1pz./h;
end