%% Lab 2: Heat Conduction Boundary Value Problem
% Solving: -d/dz(kappa * dT/dz) + v*rho*C * dT/dz = Q(z)
% Boundary conditions: T(0) = T0, T(L) = Tout

clear; clc; close all;



function [z_all, T_full] = calculateAllGridPoints(N, v)

%% ========== PARAMETERS ==========
L = 10;          % Pipe length
a = 1;           % Coil start position
b = 3;           % Coil end position
Q0 = 50;         % Heating intensity
kappa = 0.5;     % Heat conduction coefficient
rho = 1;         % Fluid density
C = 1;           % Heat capacity
T0 = 400;        % Left boundary temperature (z=0)
Tout = 300;      % Right boundary temperature (z=L)

%% ========== GRID SETUP ==========
h = L / (N + 1);              % Grid spacing
z_all = 0:h:L;                % All grid points (N+2 total, including boundaries)
z_interior = z_all(2:end-1);  % Interior points only (N points)

%% ========== COMPUTE Q(z) AT INTERIOR POINTS ==========
% Q(z) is piecewise:
%   Q = 0                           if z < a
%   Q = Q0 * sin(pi*(z-a)/(b-a))   if a <= z <= b
%   Q = 0                           if z > b

Qvec = zeros(N, 1);
for i = 1:N
    z_i = z_interior(i);
    if z_i < a
        Qvec(i) = 0;
    elseif z_i >= a && z_i <= b
        Qvec(i) = Q0 * sin(pi * (z_i - a) / (b - a));
    else
        Qvec(i) = 0;
    end


end

%% ========== DISCRETIZATION ==========
% The ODE: -kappa * d2T/dz2 + v*rho*C * dT/dz = Q(z)
%
% Using finite differences at interior point i:
%   Second derivative: d2T/dz2 ≈ (T_{i+1} - 2*T_i + T_{i-1}) / h^2
%   First derivative:  dT/dz   ≈ (T_{i+1} - T_{i-1}) / (2*h)  [central difference]
%                      or       ≈ (T_i - T_{i-1}) / h         [backward difference]
%
% Substituting into the ODE and collecting terms gives:
%   (coeff_lower) * T_{i-1} + (coeff_main) * T_i + (coeff_upper) * T_{i+1} = Q_i
%
% TODO: Work out the coefficients by hand first!
%       Substitute the finite difference formulas into the ODE,
%       then collect terms for T_{i-1}, T_i, and T_{i+1}

%% ========== BUILD TRIDIAGONAL MATRIX A ==========
% Matrix A is N x N (one row per interior point)
%
% Structure:
%   [ d   u   0   0   0  ]
%   [ l   d   u   0   0  ]
%   [ 0   l   d   u   0  ]
%   [ 0   0   l   d   u  ]
%   [ 0   0   0   l   d  ]
%
% where d = main diagonal coefficient (for T_i)
%       u = upper diagonal coefficient (for T_{i+1})
%       l = lower diagonal coefficient (for T_{i-1})

% TODO: Define your coefficients based on your discretization
coeff_lower = -(kappa/(h^2) + v*rho*C/(2*h))
coeff_main = (2*kappa/(h^2))
coeff_upper = (v*rho*C/(2*h) - kappa/(h^2))

% Build diagonal vectors
main_diag = zeros(N, 1);    % Length N
upper_diag = zeros(N-1, 1); % Length N-1 (one fewer than main)
lower_diag = zeros(N-1, 1); % Length N-1

% TODO: Fill in the diagonal vectors with your coefficients
% For uniform coefficients (same at every point):
main_diag(:) = coeff_main;
upper_diag(:) = coeff_upper;
lower_diag(:) = coeff_lower;

% Assemble the matrix using diag()
A = diag(main_diag, 0) + diag(upper_diag, 1) + diag(lower_diag, -1);

%% ========== BUILD RIGHT-HAND SIDE VECTOR b ==========
% Start with the source term
b = Qvec;

% Handle boundary conditions:
% At i=1: equation involves T_0 (known) -> move to RHS
% At i=N: equation involves T_{N+1} (known) -> move to RHS
%
% TODO: Add boundary contributions
b(1) = b(1) - (coeff_lower) * T0;
b(N) = b(N) - (coeff_upper) * Tout;

%% ========== SOLVE THE SYSTEM ==========
T_interior = A \ b;

%% ========== ASSEMBLE FULL SOLUTION ==========
% Include boundary values
T_full = [T0; T_interior; Tout];

% %% ========== PLOT RESULTS ==========
% figure;
% plot(z_all, T_full, '-o');
% xlabel('z');
% ylabel('T(z)');
% title(sprintf('Temperature Distribution (N=%d, v=%.1f)', N, v));
% grid on;
end

%% ========== CONVERGENCE STUDY (Task 1) ==========
% TODO: Run for N = 9, 19, 39, 79 with v = 0
%       Plot all curves on the same graph
%
% Hint: You can either:
%   1. Copy the code above into a loop over N values, or
%   2. Turn the code above into a function and call it multiple times
%
% Example structure:
figure;
hold on;
for N = [9, 19, 39, 79]
    [z_all, T_full] = calculateAllGridPoints(N, 0);
    plot(z_all, T_full);
end
hold off;
legend('N=9', 'N=19', 'N=39', 'N=79');

%% ========== CONVECTION STUDY (Task 2) ==========
% TODO: Use N = 79 and run for v = 0, 0.1, 0.5, 1
%       Use subplot to show all four in one figure
%
% Example structure:
figure;
v_values = [0, 0.1, 0.5, 1];
for idx = 1:4
    v = v_values(idx);
    disp(v)
    [z_all, T_full] = calculateAllGridPoints(79, v);
    subplot(2, 2, idx);
    plot(z_all, T_full);
    title(sprintf('v = %.1f', v));
    xlabel('z');
    ylabel('T(z)');
end
