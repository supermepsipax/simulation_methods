%% Lab 2: Finite Element Method - Fishing Rod Deflection
% CM2014 - Simulation Methods in Medical Engineering
%
% This script calculates the deflection of a fishing rod modeled as
% beam elements using the Finite Element Method.

clear; clc; close all;

%% ========================================================================
%  SECTION 1: Define Material and Geometric Properties
%  ========================================================================

E = 100e6; 				% Young's modulus [Pa] (100 MPa)

r = 0.05;               % Radius [m] (5 cm)
L_total = 2;            % Total rod length [m]
n_elem = 4;             % Number of elements
L_elem = L_total / n_elem;  % Length of each element [m]

A = pi * r^2;           % Cross-sectional area [m^2]
I = pi * r^4 / 4;       % Second moment of area [m^4]

beta   = 70 * pi/180;   % Angle of elements 1 and 2
lambda = 55 * pi/180;   % Angle of element 3
delta  = 30 * pi/180;   % Angle of element 4

F_applied = 20;         % Force at tip [N]

n_nodes = n_elem + 1;   % 5 nodes
n_dof = 3 * n_nodes;    % 15 DOFs (3 per node: Dx, Dz, theta)

%% ========================================================================
%  SECTION 2: Compute Node Coordinates
%  ========================================================================

nodes = zeros(n_nodes, 2);  % Column 1 = X, Column 2 = Z

nodes(1, :) = [0, 0];
nodes(2, :) = nodes(1, :) + L_elem * [cos(beta), sin(beta)];
nodes(3, :) = nodes(2, :) + L_elem * [cos(lambda), sin(lambda)];
nodes(4, :) = nodes(3, :) + L_elem * [cos(delta), sin(delta)];
nodes(5, :) = nodes(4, :) + [L_elem, 0];

elem_angles = [beta, beta, lambda, delta];

%% ========================================================================
%  SECTION 3: Helper Functions for local stiffness matrix 
%             + local -> global transformation
%  ========================================================================

% Function to compute local stiffness matrix
% Equation 7 in PDF
function ke = local_stiffness(E, A, I, L)

    ke = zeros(6, 6);

    ke(1,1) = E*A/L;
    ke(1,4) = -E*A/L;

    ke(2,2) = 12*E*I/(L^3);
    ke(2,3) = 6*E*I/(L^2);
    ke(2,5) = -12*E*I/(L^3);
    ke(2,6) = 6*E*I/(L^2);

    ke(3,2) = 6*E*I/(L^2);
    ke(3,3) = 4*E*I/(L);
    ke(3,5) = -6*E*I/(L^2);
    ke(3,6) = 2*E*I/(L);

    ke(4,1) = -E*A/L;
    ke(4,4) = E*A/L;

    ke(5,2) = -12*E*I/(L^3);
    ke(5,3) = -6*E*I/(L^2);
    ke(5,5) = 12*E*I/(L^3);
    ke(5,6) = -6*E*I/(L^2);

    ke(6,2) = 6*E*I/(L^2);
    ke(6,3) = 2*E*I/(L);
    ke(6,5) = -6*E*I/(L^2);
    ke(6,6) = 4*E*I/(L);
end

function T = transformation_matrix(angle)
    c = cos(angle);
	s = sin(angle);

    T = zeros(6, 6);

    T(1,1) = c;
	T(1,2) = s;   
    
	T(2,1) = -s;  
	T(2,2) = c;   

	T(3,3) = 1;

    T(4,4) = c;
	T(4,5) = s;

    T(5,4) = -s;
	T(5,5) = c;
    
	T(6,6) = 1;
end

%% ========================================================================
%  SECTION 4: Assemble Global Stiffness Matrix
%  ========================================================================

% global stiffness matrix
K = zeros(n_dof, n_dof);

for e = 1:n_elem
    start_node = e;
    end_node = e + 1;

    angle = elem_angles(e);

    ke = local_stiffness(E, A, I, L_elem);

    T = transformation_matrix(angle);

    Ke = T' * ke * T;

    % global DOF indices
    dof = [3*(start_node-1)+(1:3), 3*(end_node-1)+(1:3)];

    % add element matrix to global matrix
    K(dof, dof) = K(dof, dof) + Ke;
end

%% ========================================================================
%  SECTION 5: Apply Boundary Conditions and Solve
%  ========================================================================

% Create force vector
F = zeros(n_dof, 1);

F(14) = -F_applied;

% Apply boundary conditions (Node 1 is fixed: D1 = D2 = D3 = 0)
fixed_dof = [1, 2, 3];

free_dof = setdiff(1:n_dof, fixed_dof);

K_reduced = K(free_dof, free_dof);
F_reduced = F(free_dof);
D_free = K_reduced \ F_reduced;

D = zeros(n_dof, 1);
D(free_dof) = D_free;

%% ========================================================================
%  SECTION 6: Post-Processing and Results
%  ========================================================================

% Reshape displacement vector for easier access
% Each row = one node, columns = [Dx, Dz, theta]
D_reshaped = reshape(D, 3, n_nodes)';

Dx = D_reshaped(:, 1);      % X-displacements
Dz = D_reshaped(:, 2);      % Z-displacements
theta = D_reshaped(:, 3);   % Rotations

X_def = nodes(:, 1) + Dx;
Z_def = nodes(:, 2) + Dz;

D_total = sqrt(Dx.^2 + Dz.^2);

%% ========================================================================
%  SECTION 7: Display Results Table
%  ========================================================================

fprintf('\n========================================\n');
fprintf('   NODAL DISPLACEMENTS\n');
fprintf('========================================\n');
fprintf('Node\t Dx (m)\t\t Dz (m)\t\t Theta (rad)\t Total (m)\n');
fprintf('----\t-------\t\t-------\t\t-----------\t---------\n');

for i = 1:n_nodes
    fprintf('%d\t%10.6f\t%10.6f\t%10.6f\t%10.6f\n', ...
            i, Dx(i), Dz(i), theta(i), D_total(i));
end

fprintf('\nMaximum total displacement: %.4f m (at node %d)\n', ...
        max(D_total), find(D_total == max(D_total)));

%% ========================================================================
%  SECTION 8: Visualization
%  ========================================================================

figure('Position', [100, 100, 800, 600]);

% --- Plot undeformed and deformed shapes ---
subplot(1, 1, 1);
hold on;

% Undeformed shape
plot(nodes(:,1), nodes(:,2), 'b-o', 'LineWidth', 2, 'MarkerSize', 10, ...
     'MarkerFaceColor', 'b', 'DisplayName', 'Undeformed');

% Deformed shape
% Adjust scale factor if displacements are too small to see clearly
scale = 1;  % Set to 1 for actual displacement, increase to exaggerate
plot(X_def, Z_def, 'r-s', 'LineWidth', 2, 'MarkerSize', 10, ...
     'MarkerFaceColor', 'r', 'DisplayName', sprintf('Deformed (scale=%.0f)', scale));

% Mark the fixed support
plot(0, 0, 'k^', 'MarkerSize', 15, 'MarkerFaceColor', 'k');

% Draw the force arrow at the tip
% TODO: Add an arrow showing the applied force direction (optional)

% Formatting
legend('Location', 'best');
xlabel('X (m)');
ylabel('Z (m)');
title('Fishing Rod Deflection Under 20 N Load');
axis equal;
grid on;

% Adjust axis limits for better visibility
xlim([-0.1, max(nodes(:,1)) + 0.15]);
ylim([-0.1, max(nodes(:,2)) + 0.15]);

hold off;

%% ========================================================================
%  SECTION 9: Verification (Optional)
%  ========================================================================


fprintf('\n========================================\n');
fprintf('   VERIFICATION\n');
fprintf('========================================\n');
fprintf('Expected (from Comsol):\n');
fprintf('  Max total displacement: ~0.0796 m\n');
fprintf('  Tip X-displacement: ~0.0361 m\n');
fprintf('  Tip Z-displacement: ~0.071 m\n');
fprintf('\nYour results:\n');
fprintf('  Max total displacement: %.4f m\n', max(D_total));
fprintf('  Tip X-displacement: %.4f m\n', Dx(end));
fprintf('  Tip Z-displacement: %.4f m\n', Dz(end));
