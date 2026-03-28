
% 历时 78.343651 秒。
% 历时 70.270201 秒。
clear; clc; close all

% --- Parameters ---
nt = 2000;
nx = 450; nz = 450;
dx = 15; h = dx;
dt = 0.0021; 

% Velocity model
v = ones(nz,nx)*3500;
v(1:nz/2,:) = 2200;

% Source
f0 = 20*pi;
t = (0:nt-1)*dt;
t0 = 4/f0;
src_raw = 10^6 * exp(-f0^2*(t-t0).^2);
src = -[0, diff(src_raw)]/dx^2;

% Fields
P = zeros(nz,nx); Vx = P; Vz = P;
zs = 51; xs = round(nx/2);

% High-order coefficients (M=7)
coeff = [1.55427, -0.27162, 0.07394, -0.0208847, 0.005117, -0.000919339, 0.0000882575];
ng = 7; % ghost width

% RKS coefficients
% d = [22/24, 1/24, 1/24];
d = [0.833333, 0.0833333, 0.0833333];

d_prime = d;

seis_record = zeros(nt,nx);
tic
for it = 1:nt-2
  
    [dvx_dx, dvz_dz] = div_M7_FS(Vx, Vz, coeff, h, ng);
    P = P + dt * v.^2 .* (dvx_dx + dvz_dz) ;

    P(zs,xs) = P(zs,xs) + src(it);
    P(1,:) = 0; % Free Surface at top row
    P = spongeABC_noTop(P, nx, nz, 50, 0.009);


%%%%%%%%%%%%%%%%%%


    p2_prime = P - dt *v.^2.* (dvx_dx + dvz_dz);
 
       [dpx_dx, dpy_dy] = grad_M1_FS(P,h);
    vx2_prime = Vx + dt * dpx_dx;
    vy2_prime = Vz + dt * dpy_dy;

      [dvx_dx, dvy_dy] = div_M7_FS(vx2_prime, vy2_prime,coeff,h,ng);
    p3_prime = P + dt *v.^2.* (dvx_dx + dvy_dy);

       [dpx_dx1, dpy_dy1] = grad_M1_FS(p2_prime,h);
    [dpx_dx2, dpy_dy2] = grad_M1_FS(p3_prime,h);

        Vx = Vx + dt * (d_prime(1)*dpx_dx + ...
        d_prime(2)*dpx_dx1 + ...
        d_prime(3)*dpx_dx2);

    Vz = Vz + dt * (d_prime(1)*dpy_dy + ...
        d_prime(2)*dpy_dy1 + ...
        d_prime(3)*dpy_dy2);


    Vx = spongeABC_noTop(Vx, nx, nz, 50, 0.009);
    Vz = spongeABC_noTop(Vz, nx, nz, 50, 0.009);

    if rem(it,60) == 0
        imagesc(P, [-1 1]*50); colormap gray; axis equal;
        title(sprintf('RKS Free Surface | Step: %i', it));
        drawnow;
    end
    seis_record(it,:) = P(zs,:);
end
toc
save('figure2d.mat')
%% --- FUNCTIONS ---

% function [dpx, dpz] = grad_M1_FS(P, h)
%     % 2nd order staggered gradient (Pressure to Velocity)
%     [nz, nx] = size(P);
%     dpx = zeros(nz, nx); dpz = zeros(nz, nx);
%     % dP/dx: Forward difference
%     dpx(:, 1:nx-1) = (P(:, 2:nx) - P(:, 1:nx-1))/h;
%     % dP/dz: Forward difference. P(1)=0 is surface, so dpz(1) is grad at depth h/2
%     dpz(1:nz-1, :) = (P(2:nz, :) - P(1:nz-1, :))/h;
% end

function [dpx,dpz] =grad_M1_FS(p, h)

dpx=(circshift(p,[0 -1])-circshift(p,[0 0]));
dpz=(circshift(p,[-1])-circshift(p,[0]));
dpx=dpx./h;
dpz=dpz./h;
end

function [dvx_dx, dvz_dz] = div_M7_FS(Vx, Vz, coeff, h, ng)
    % 7th order staggered divergence (Velocity to Pressure)
    [nz, nx] = size(Vx);
    dvx_dx = zeros(nz, nx); dvz_dz = zeros(nz, nx);
    
    % Padding for high-order reach
    % Vz uses EVEN symmetry at top for Free Surface
    Vzg = padarray(Vz, [ng ng], 'replicate');
    for k = 1:ng
        Vzg(ng+1-k, :) = Vzg(ng+1+k, :); 
    end
    Vxg = padarray(Vx, [ng ng], 'replicate');

    iz = ng+1:ng+nz; ix = ng+1:ng+nx;
    for m = 1:7
        dvx_dx = dvx_dx + coeff(m) * (Vxg(iz, ix+m-1) - Vxg(iz, ix-m));
        dvz_dz = dvz_dz + coeff(m) * (Vzg(iz+m-1, ix) - Vzg(iz-m, ix));
    end
    dvx_dx = dvx_dx / h; dvz_dz = dvz_dz / h;
end

function field = spongeABC_noTop(field, nx, nz, span, alpha)
    % Damping on sides and bottom only to allow top reflection
    c2 = alpha^2;
    for t = 1:span
        w = exp(-c2*(span-t)^2);
        field(:, t) = field(:, t) * w;            % Left
        field(:, nx-t+1) = field(:, nx-t+1) * w;  % Right
        field(nz-t+1, :) = field(nz-t+1, :) * w;  % Bottom
    end
end