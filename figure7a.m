% 历时 653.102796 秒。
% 历时 652.022056 秒。
clear
clc %%%%%%%
close all
% Elapsed time is 9.253353 seconds.
nt=3539;    % number of time steps
eps=.6;     % stability
isnap=60;    % snapshot sampling

load('vv.mat')
c1=flipud(c);
c=c1;

[nz,nx]=size(c1);
v=c1;


nz=nz+50;

vv=zeros(nz,nx);
for ii=1:nz-50
    for jj=1:nx
        vv(ii+50,jj)=v(ii,jj);
    end
end

for ii=1:50  %%top
    for jj=1:nx
        vv(ii,jj)=vv(51,jj);
    end
end



clear v
v=vv;



p=zeros([nz nx]); pboundarynew=p;pdan=p;
dx=15;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.0015; % calculate time step from stability criterion
tau=dt;
r2=v.*v.*dt*dt/dx/dx;


f0=45;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^6*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain

seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=p;
d2pz=p;

% Source location
% Source location
zs=51;
xs=nx/2;

h=dx;
r=v*dt/h;

p=zeros([nz nx]); Vx=p; Vz=p;
coeff=[ 1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];

d = [22/24, 1/24, 1/24];
d_prime=d;

dpx_dx=zeros(nz,nx);
dpy_dy=zeros(nz,nx);
Vx2=zeros(nz,nx);
Vz2=zeros(nz,nx);
tic
for it=1:nt-2,

    [dvx_dx, dvy_dy] = Vx_Vz_spatial_derivatives(Vx,Vz,coeff);
    p2=p+dt*v.^2.*(dvx_dx+dvy_dy);

    [dpx_dx2, dpy_dy2] = spatial_derivatives(p2,coeff);
    Vx3 = Vx + dt * dpx_dx2;
    Vz3 = Vz + dt * dpy_dy2;


    % Stage 5
    Vx2 = Vx - dt * dpx_dx;
    Vz2 = Vz - dt * dpy_dy;
    [dvx_dx1, dvy_dy1] = Vx_Vz_spatial_derivatives(Vx2, Vz2,coeff);
    [dvx_dx2, dvy_dy2] = Vx_Vz_spatial_derivatives(Vx3, Vz3,coeff);


    % p=p-dt*v.^2.*(Vxd1+Vzd1)/h;
    p = p + dt * v.^2.*(d(1).*(dvx_dx + dvy_dy) + ...
        d(2).*(dvx_dx1 + dvy_dy1) + ...
        d(3).*(dvx_dx2 + dvy_dy2));

    p(zs,xs)= p(zs,xs)+src(it);
    [p,p]=spongeABC(p,p,nx,nz,50,50,0.009);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    p2_prime = p - dt *v.^2.* (dvx_dx + dvy_dy);

    % Stage 3'
    [dpx_dx, dpy_dy] = spatial_derivatives(p,coeff);
    vx2_prime = Vx + dt * dpx_dx;
    vy2_prime = Vz + dt * dpy_dy;

    % Stage 4'
    [dvx_dx, dvy_dy] = Vx_Vz_spatial_derivatives(vx2_prime, vy2_prime,coeff);
    % [dvx_dx, dvy_dy] = Vx_Vz_spatial_derivatives_short(vx2_prime, vy2_prime);
    p3_prime = p + dt *v.^2.* (dvx_dx + dvy_dy);

    % Stage 5'
    [dpx_dx1, dpy_dy1] = spatial_derivatives(p2_prime,coeff);
    [dpx_dx2, dpy_dy2] = spatial_derivatives(p3_prime,coeff);

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
save('figure7a_TraRK.mat')
function [d1px,d1pz] =spatial_derivatives(p,coeff)
h=15;
d1px1=(circshift(p,[0 -1])-circshift(p,[0 0]));
d1px2=(circshift(p,[0 -2])-circshift(p,[0 1]));
d1px3=(circshift(p,[0 -3])-circshift(p,[0 2]));
d1px4=(circshift(p,[0 -4])-circshift(p,[0 3]));
d1px5=(circshift(p,[0 -5])-circshift(p,[0 4]));
d1px6=(circshift(p,[0 -6])-circshift(p,[0 5]));
d1px7=(circshift(p,[0 -7])-circshift(p,[0 6]));


d1pz1=(circshift(p,[-1])-circshift(p,[0]));
d1pz2=(circshift(p,[-2])-circshift(p,[1]));
d1pz3=(circshift(p,[-3])-circshift(p,[2]));
d1pz4=(circshift(p,[-4])-circshift(p,[3]));
d1pz5=(circshift(p,[-5])-circshift(p,[4]));
d1pz6=(circshift(p,[-6])-circshift(p,[5]));
d1pz7=(circshift(p,[-7])-circshift(p,[6]));

d1px=coeff(1)*d1px1+coeff(2)*d1px2+coeff(3)*d1px3+coeff(4)*d1px4+coeff(5)*d1px5+coeff(6)*d1px6...
    +coeff(7)*d1px7;
d1pz=coeff(1)*d1pz1+coeff(2)*d1pz2+coeff(3)*d1pz3+coeff(4)*d1pz4+coeff(5)*d1pz5+coeff(6)*d1pz6...
    +coeff(7)*d1pz7;
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