% . (Modified for Free Surface Condition)

% 历时 79.042632 秒。
clear; clc; close all

% --- Parameters ---
nt = floor(2000*2.1/6.1);
nx = 450; nz = 450;
dx = 15; h = dx;
dt = 0.0061;

% Velocity model
v = ones(nz,nx)*3500;
v(1:nz/2,:) = 2200;

% Source
f0 = 20*pi;
t = (0:nt-1)*dt;
t0 = 4/f0;
src_func = 10^6 * exp(-f0^2*(t-t0).^2);
src = -[0, diff(src_func)]/dx^2;

% Fields
P = zeros(nz,nx); Vx = P; Vz = P;
zs = 51; xs = round(nx/2);

% High-order coefficients (M=7)
coeff = [1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];
ng = 7; % ghost width

% RKS coefficients
d = [22/24, 1/24, 1/24];
d_prime = d;

seis_record = zeros(nt,nx);

dpx_dx_p=P;
dpz_dz_p=P;
tic
for it = 1:nt-2
    % --- STAGE 1: PRESSURE UPDATES ---
    [dvx_dx, dvz_dz] = div_M7_FS(Vx, Vz, coeff, h, ng);

    % Predictor P2
    p2 = P + dt * v.^2 .* (dvx_dx + dvz_dz);
    p2(1,:) = 0; % Enforce Free Surface

    [dpx_dx2, dpz_dz2] = grad_M7_FS(p2, coeff, h, ng);
    Vx3 = Vx + dt * dpx_dx2;
    Vz3 = Vz + dt * dpz_dz2;

    % [dpx_now, dpz_now] = grad_M7_FS(P, coeff, h, ng);
    Vx2 = Vx - dt * dpx_dx_p;
    Vz2 = Vz - dt * dpz_dz_p;

    [dvx_dx1, dvz_dz1] = div_M7_FS(Vx2, Vz2, coeff, h, ng);
    [dvx_dx2, dvz_dz2] = div_M7_FS(Vx3, Vz3, coeff, h, ng);

    % Final Pressure Update (RKS Combination)
    P = P + dt * v.^2 .* (d(1)*(dvx_dx + dvz_dz) + ...
        d(2)*(dvx_dx1 + dvz_dz1) + ...
        d(3)*(dvx_dx2 + dvz_dz2));

    P(zs,xs) = P(zs,xs) + src(it);
    P(1,:) = 0; % Ensure Pressure release at top boundary

    % Sponge ABC (Sides and Bottom only)
    P = spongeABC_noTop(P, nx, nz, 50, 0.009);

    % --- STAGE 2: VELOCITY UPDATES ---
    p2_prime = P - dt * v.^2 .* (dvx_dx + dvz_dz);
    p2_prime(1,:) = 0;

    [dpx_dx_p, dpz_dz_p] = grad_M7_FS(P, coeff, h, ng);
    vx2_prime = Vx + dt * dpx_dx_p;
    vz2_prime = Vz + dt * dpz_dz_p;

    [dvx_dx_p3, dvz_dz_p3] = div_M7_FS(vx2_prime, vz2_prime, coeff, h, ng);
    p3_prime = P + dt * v.^2 .* (dvx_dx_p3 + dvz_dz_p3);
    p3_prime(1,:) = 0;

    [dpx_dx_s1, dpz_dz_s1] = grad_M7_FS(p2_prime, coeff, h, ng);
    [dpx_dx_s2, dpz_dz_s2] = grad_M7_FS(p3_prime, coeff, h, ng);

    Vx = Vx + dt * (d_prime(1)*dpx_dx_p  + d_prime(2)*dpx_dx_s1  + d_prime(3)*dpx_dx_s2);
    Vz = Vz + dt * (d_prime(1)*dpz_dz_p  + d_prime(2)*dpz_dz_s1  + d_prime(3)*dpz_dz_s2);

    Vx = spongeABC_noTop(Vx, nx, nz, 50, 0.009);
    Vz = spongeABC_noTop(Vz, nx, nz, 50, 0.009);

    if rem(it,60) == 0
        imagesc(P, [-1 1]*50); colormap gray; axis equal;
        title(sprintf('High-Order RKS Free Surface | Step: %i', it));
        drawnow;
    end
    seis_record(it,:) = P(zs,:);
end
toc
save('figure5a.mat');

%% --- FUNCTIONS ---

function [dpx, dpz] = grad_M7_FS(P, coeff, h, ng)
% Staggered gradient (Pressure to Velocity) with Odd Symmetry top
[nz, nx] = size(P);
dpx = zeros(nz, nx); dpz = zeros(nz, nx);

% Pad Pressure: ODD symmetry at top ( ghost = -P ), replicate others
Pg = padarray(P, [ng ng], 'replicate');
for k = 1:ng
    Pg(ng+1-k, :) = -Pg(ng+1+k, :); % Odd Mirror for Free Surface
end

iz = ng+1:ng+nz; ix = ng+1:ng+nx;
for m = 1:7
    dpx = dpx + coeff(m)*(Pg(iz, ix+m) - Pg(iz, ix-m+1));
    dpz = dpz + coeff(m)*(Pg(iz+m, ix) - Pg(iz-m+1, ix));
end
dpx = dpx/h; dpz = dpz/h;
end

function [dvx_dx, dvz_dz] = div_M7_FS(Vx, Vz, coeff, h, ng)
% Staggered divergence (Velocity to Pressure) with Even Symmetry top
[nz, nx] = size(Vx);
dvx_dx = zeros(nz, nx); dvz_dz = zeros(nz, nx);

% Pad Velocity: Vz needs EVEN symmetry at top, replicate others
Vzg = padarray(Vz, [ng ng], 'replicate');
for k = 1:ng
    Vzg(ng+1-k, :) = Vzg(ng+1+k, :); % Even Mirror
end
Vxg = padarray(Vx, [ng ng], 'replicate');

iz = ng+1:ng+nz; ix = ng+1:ng+nx;
for m = 1:7
    dvx_dx = dvx_dx + coeff(m)*(Vxg(iz, ix+m-1) - Vxg(iz, ix-m));
    dvz_dz = dvz_dz + coeff(m)*(Vzg(iz+m-1, ix) - Vzg(iz-m, ix));
end
dvx_dx = dvx_dx/h; dvz_dz = dvz_dz/h;
end

function field = spongeABC_noTop(field, nx, nz, span, alpha)
% Damping on Left, Right, and Bottom. Top is skipped for Free Surface.
c2 = alpha^2;
for t = 1:span
    w = exp(-c2*(span-t)^2);
    field(:, t) = field(:, t) * w;            % Left
    field(:, nx-t+1) = field(:, nx-t+1) * w;  % Right
    field(nz-t+1, :) = field(nz-t+1, :) * w;  % Bottom
end
end