% 历时 62.774831 秒。 (TDE + High-Order SGFD + Free Surface)
clear; clc; close all

% --- 1. Basic Parameters ---
nt = 2000;
nx = 450; nz = 450;
dx = 15; h = dx;
dt = 0.0021; 

% Velocity model
v = ones(nz,nx)*3500;
v(1:nz/2,:) = 2200;

% Grid coordinates
x = (0:nx-1)*dx;
z = (0:nz-1)*dx;

% --- 2. Source with TDE (Temporal Dispersion Transform) ---
f0 = 20*pi;
t_axis = (0:nt-1)*dt;
t0 = 4/f0;
src_raw = 10^6 * exp(-f0^2 * (t_axis - t0).^2);
src_raw = -[0, diff(src_raw)]/dx^2; % First derivative

% TDE Matrix Transformation
nt_tde = nt - 1;
B = ifft(exp(-2i*sin([0:nt_tde-1]*pi/(2*nt_tde))'*[0:nt_tde-1]), 2*nt_tde, 'symmetric');
T = B(1:nt_tde, 1:nt_tde); 
src = T * src_raw(1:nt_tde)'; % Transforming source

% --- 3. Field Initialization ---
P = zeros(nz,nx); Vx = P; Vz = P;
zs = 51; xs = round(nx/2); % Source position

% High-order coefficients (M=7)
coeff = [1.24275, -0.113493, 0.0272555, -0.00731429, 0.00173544, -0.000303651, 0.0000284348];

seis_record = zeros(nt_tde, nx);
tic
for it = 1:nt_tde - 1
    
    % --- Step A: Update Pressure (Divergence of V) ---
    [dvx_dx, dvz_dz] = div_M7_FS(Vx, Vz, coeff, h);
    P = P + dt * (v.^2) .* (dvx_dx + dvz_dz);

    % Inject TDE-transformed source
    P(zs, xs) = P(zs, xs) + src(it);

    % FREE SURFACE: Explicitly zero the top row
    P(1,:) = 0;

    % Sponge ABC (Sides and Bottom ONLY)
    P = spongeABC_noTop(P, nx, nz, 50, 0.009);

    % --- Step B: Update Velocity (Gradient of P) ---
    [dpx_dx, dpz_dz] = grad_M7_FS(P, coeff, h);
    Vx = Vx + dt * dpx_dx;
    Vz = Vz + dt * dpz_dz;

    % Sponge ABC for Velocity
    Vx = spongeABC_noTop(Vx, nx, nz, 50, 0.009);
    Vz = spongeABC_noTop(Vz, nx, nz, 50, 0.009);

    % --- Visualization ---
    if rem(it, 60) == 0
        imagesc(x, z, P, [-1 1]*50); axis equal;
        colormap gray; colorbar;
        title(sprintf('TDE + Free Surface | Step: %i', it));
        xlabel('x (m)'); ylabel('z (m)');
        drawnow;
    end
    seis_record(it,:) = P(zs,:);
end
toc

% save('figure2f_TDE_FreeSurface.mat', 'seis_record');
save('figure2f.mat')

%% --- FUNCTIONS ---

function [dvx_dx, dvz_dz] = div_M7_FS(Vx, Vz, coeff, h)
    [nz, nx] = size(Vx);
    dvx_dx = zeros(nz, nx); dvz_dz = zeros(nz, nx);
    
    % Padding for high-order reach (7 cells)
    % Vz uses EVEN symmetry at the top (z=0)
    Vzg = padarray(Vz, [7 7], 'replicate');
    for k = 1:7
        Vzg(8-k, :) = Vzg(8+k, :); % Mirror top
    end
    Vxg = padarray(Vx, [7 7], 'replicate'); % Standard replicate for sides/bottom

    iz = 8:8+nz-1; ix = 8:8+nx-1;
    for m = 1:7
        % Divergence Staggering Logic
        dvx_dx = dvx_dx + coeff(m) * (Vxg(iz, ix+m-1) - Vxg(iz, ix-m));
        dvz_dz = dvz_dz + coeff(m) * (Vzg(iz+m-1, ix) - Vzg(iz-m, ix));
    end
    dvx_dx = dvx_dx/h; dvz_dz = dvz_dz/h;
end

function [dpx_dx, dpz_dz] = grad_M7_FS(P, coeff, h)
    [nz, nx] = size(P);
    dpx_dx = zeros(nz, nx); dpz_dz = zeros(nz, nx);
    
    % Padding for high-order reach
    % Pressure uses ODD symmetry at the top (P=0)
    Pg = padarray(P, [7 7], 'replicate');
    for k = 1:7
        Pg(8-k, :) = -Pg(8+k, :); % Mirror top with negative polarity
    end
    
    iz = 8:8+nz-1; ix = 8:8+nx-1;
    for m = 1:7
        % Gradient Staggering Logic
        dpx_dx = dpx_dx + coeff(m) * (Pg(iz, ix+m) - Pg(iz, ix-m+1));
        dpz_dz = dpz_dz + coeff(m) * (Pg(iz+m, ix) - Pg(iz-m+1, ix));
    end
    dpx_dx = dpx_dx/h; dpz_dz = dpz_dz/h;
end

function field = spongeABC_noTop(field, nx, nz, span, alpha)
    % Damping restricted to Left, Right, and Bottom boundaries
    % Top boundary (z=1 to span) is ignored to preserve reflection
    c2 = alpha^2;
    % Left/Right
    for t = 1:span
        w = exp(-c2*(span-t)^2);
        field(:, t) = field(:, t) * w;            % Left
        field(:, nx-t+1) = field(:, nx-t+1) * w;  % Right
    end
    % Bottom
    for t = nz-span+1:nz
        w = exp(-c2*(t-nz+span-1)^2);
        field(t, :) = field(t, :) * w;
    end
end