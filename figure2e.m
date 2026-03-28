% Modified High-Order SGFD with Free Surface Condition
% Based on Liang & Wang Manuscript Logic
clear; clc; close all

% --- Parameters ---
nt = 2000;
nx = 450; nz = 450;
dx = 15; h = dx;
dt = 0.0021; % CFL check: v_max*dt/dx = 3500*0.0021/15 = 0.49 (Stable)

% Velocity model
v = ones(nz,nx)*3500;
v(1:nz/2,:) = 2200;

% Source (Ricker derivative)
f0 = 20*pi;
t = (0:nt-1)*dt;
t0 = 4/f0;
src_func = 10^6 * exp(-f0^2*(t-t0).^2);
src = -[0, diff(src_func)]/dx^2;

% Source location
zs = 51; xs = round(nx/2);

% Fields
P = zeros(nz,nx); Vx = P; Vz = P;
seis_record = zeros(nt,nx);

% High-order coefficients (M=7)
coeff = [1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];

tic
for it = 1:nt-2
    % --- 1. Update Pressure (Using Divergence of V) ---
    [dvx_dx, dvz_dz] = div_M7_FreeSurface(Vx, Vz, coeff, h);
    P = P + dt * (v.^2) .* (dvx_dx + dvz_dz);

    % Inject Source
    P(zs, xs) = P(zs, xs) + src(it);

    % --- FREE SURFACE CONDITION ---
    P(1,:) = 0; % Pressure release at top row

    % Sponge ABC (Modified: No damping at top row to allow reflection)
    P = spongeABC_noTop(P, nx, nz, 50, 0.009);

    % --- 2. Update Velocities (Using Gradient of P) ---
    [dpx_dx, dpz_dz] = grad_M7_FreeSurface(P, coeff, h);
    Vx = Vx + dt * dpx_dx;
    Vz = Vz + dt * dpz_dz;

    % Sponge ABC for Velocity
    Vx = spongeABC_noTop(Vx, nx, nz, 50, 0.009);
    Vz = spongeABC_noTop(Vz, nx, nz, 50, 0.009);

    % --- Visualization ---
    if rem(it, 60) == 0
        imagesc(P, [-1 1]*50); axis equal; colormap gray;
        title(sprintf('Step: %i | Free Surface (P=0) at Top', it));
        drawnow;
    end
    seis_record(it,:) = P(zs,:);
end
toc

save('figure2e_2order_SGFD.mat')

%% --- HELPER FUNCTIONS FOR FREE SURFACE ---

function [dvx_dx, dvz_dz] = div_M7_FreeSurface(Vx, Vz, coeff, h)
    [nz, nx] = size(Vx);
    dvx_dx = zeros(nz, nx);
    dvz_dz = zeros(nz, nx);
    
    % Padding to handle M=7 reach without circshift wrap
    % Vz uses Even Symmetry for the top boundary
    Vzg = padarray(Vz, [7 7], 'replicate');
    % Mirror top for Vz (Even symmetry: Vz(-z) = Vz(z))
    for k = 1:7
        Vzg(8-k, :) = Vzg(8+k, :); 
    end
    Vxg = padarray(Vx, [7 7], 'replicate');

    iz = 8:8+nz-1; ix = 8:8+nx-1;
    for m = 1:7
        dvx_dx = dvx_dx + coeff(m) * (Vxg(iz, ix+m-1) - Vxg(iz, ix-m));
        dvz_dz = dvz_dz + coeff(m) * (Vzg(iz+m-1, ix) - Vzg(iz-m, ix));
    end
    dvx_dx = dvx_dx / h;
    dvz_dz = dvz_dz / h;
end

function [dpx_dx, dpz_dz] = grad_M7_FreeSurface(P, coeff, h)
    [nz, nx] = size(P);
    dpx_dx = zeros(nz, nx);
    dpz_dz = zeros(nz, nx);
    
    % Pressure uses Odd Symmetry for the top boundary (P=0)
    Pg = padarray(P, [7 7], 'replicate');
    % Mirror top for P (Odd symmetry: P(-z) = -P(z))
    for k = 1:7
        Pg(8-k, :) = -Pg(8+k, :); 
    end
    
    iz = 8:8+nz-1; ix = 8:8+nx-1;
    for m = 1:7
        dpx_dx = dpx_dx + coeff(m) * (Pg(iz, ix+m) - Pg(iz, ix-m+1));
        dpz_dz = dpz_dz + coeff(m) * (Pg(iz+m, ix) - Pg(iz-m+1, ix));
    end
    dpx_dx = dpx_dx / h;
    dpz_dz = dpz_dz / h;
end

function field = spongeABC_noTop(field, nx, nz, span, coeff_damp)
    % Standard exponential sponge, but skipping the top boundary (1:span in Z)
    % to prevent dampening the free surface reflection.
    coeff2 = coeff_damp^2;
    % Left
    for t = 1:span
        w = exp(-coeff2*(span-t)^2);
        field(:, t) = field(:, t) * w;
    end
    % Right
    for t = nx-span+1:nx
        w = exp(-coeff2*(t-nx+span-1)^2);
        field(:, t) = field(:, t) * w;
    end
    % Bottom
    for t = nz-span+1:nz
        w = exp(-coeff2*(t-nz+span-1)^2);
        field(t, :) = field(t, :) * w;
    end
end