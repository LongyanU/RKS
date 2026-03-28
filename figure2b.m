
% 历时 144.613051 秒。
% 历时 144.905929 秒。
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
src_func = 10^6 * exp(-f0^2*(t-t0).^2);
src = -[0, diff(src_func)]/dx^2;

% Fields
P = zeros(nz,nx); Vx = P; Vz = P;
zs = 51; xs = round(nx/2);

% High-order coefficients (M=7)
coeff = [1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];
ng = 7; % ghost width

% RKS coefficients for NU-T scheme (Equation 19)
d = [0.833333, 0.0833333, 0.0833333];
dpz_final=P;
dpx_final=P;
seis_record = zeros(nt,nx);
tic
for it = 1:nt-2
    % --- STAGE 1: PRESSURE UPDATES (NU-T Logic) ---
    % High-order divergence of current velocities
    [dvx_dx, dvz_dz] = div_M7_FS(Vx, Vz, coeff, h, ng);
    
    % Predictor P2
    p2 = P + dt * v.^2 .* (dvx_dx + dvz_dz);
    p2(1,:) = 0; % Enforce Free Surface on intermediate field
    
    % Gradients of P2 
    [dpx_dx2, dpz_dz2] = grad_M7_FS(p2, coeff, h, ng);
    Vx3 = Vx + dt * dpx_dx2;
    Vz3 = Vz + dt * dpz_dz2;

    % Intermediate stages

    Vx2 = Vx - dt * dpx_final;
    Vz2 = Vz - dt * dpz_final;
    
    [dvx_dx1, dvz_dz1] = div_M7_FS(Vx2, Vz2, coeff, h, ng);
    [dvx_dx2, dvz_dz2] = div_M7_FS(Vx3, Vz3, coeff, h, ng);

    % Final Pressure Update (RKS Combination)
    P = P + dt * v.^2 .* (d(1)*(dvx_dx + dvz_dz) + ...
                          d(2)*(dvx_dx1 + dvz_dz1) + ...
                          d(3)*(dvx_dx2 + dvz_dz2));

    % Inject Source and Finalize P Boundary
    P(zs,xs) = P(zs,xs) + src(it);
    P(1,:) = 0; % Force Pressure Release at top row
    
    % Sponge ABC (Sides and Bottom only)
    P = spongeABC_noTop(P, nx, nz, 50, 0.009);

    % --- STAGE 2: VELOCITY UPDATE (Standard Update for NU-T) ---
    [dpx_final, dpz_final] = grad_M7_FS(P, coeff, h, ng);
    Vx = Vx + dt * dpx_final;
    Vz = Vz + dt * dpz_final;

    % Sponge ABC for Velocity
    Vx = spongeABC_noTop(Vx, nx, nz, 50, 0.009);
    Vz = spongeABC_noTop(Vz, nx, nz, 50, 0.009);

    % --- Visualization ---
    if rem(it, 60) == 0
        imagesc(P, [-1 1]*50); colormap gray; axis equal;
        title(sprintf('NU-T-RKS Free Surface | Step: %i', it));
        xlabel('X'); ylabel('Z');
        drawnow;
    end
    seis_record(it,:) = P(zs,:);
end
toc
save('figure2b.mat');

%% --- FUNCTIONS ---

function [dpx, dpz] = grad_M7_FS(P, coeff, h, ng)
    % High-order staggered gradient (Pressure to Velocity)
    [nz, nx] = size(P);
    dpx = zeros(nz, nx); dpz = zeros(nz, nx);
    
    % Padding Pressure with ODD symmetry at top (forces P=0 at boundary)
    Pg = padarray(P, [ng ng], 'replicate');
    for k = 1:ng
        Pg(ng+1-k, :) = -Pg(ng+1+k, :); % Odd Mirror
    end
    
    iz = ng+1:ng+nz; ix = ng+1:ng+nx;
    for m = 1:7
        dpx = dpx + coeff(m)*(Pg(iz, ix+m) - Pg(iz, ix-m+1));
        dpz = dpz + coeff(m)*(Pg(iz+m, ix) - Pg(iz-m+1, ix));
    end
    dpx = dpx/h; dpz = dpz/h;
end

function [dvx_dx, dvz_dz] = div_M7_FS(Vx, Vz, coeff, h, ng)
    % High-order staggered divergence (Velocity to Pressure)
    [nz, nx] = size(Vx);
    dvx_dx = zeros(nz, nx); dvz_dz = zeros(nz, nx);
    
    % Padding Velocity: Vz needs EVEN symmetry at the top
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
    % Damping on Left, Right, and Bottom. Top is excluded for Free Surface.
    c2 = alpha^2;
    for t = 1:span
        w = exp(-c2*(span-t)^2);
        field(:, t) = field(:, t) * w;            % Left
        field(:, nx-t+1) = field(:, nx-t+1) * w;  % Right
        field(nz-t+1, :) = field(nz-t+1, :) * w;  % Bottom
    end
end