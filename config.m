function [fd, params, dof] = config()
%CONFIG Build model parameters, mesh data, and finite-difference operators.
%   [fd, params, dof] = config() returns the geometry, solver parameters,
%   and discretization metadata used by the steady-state solver, the
%   concentration-coupled solver, and the visualization scripts.
%
%   Outputs:
%     fd     Mesh coordinates and finite-difference operators.
%     params Model, runtime, and concentration parameters.
%     dof    Discretization sizes and save cadence metadata.
%
%   Notes:
%     - The current solver assumes the cell ordering
%       [nurse cells, oocyte, cluster].
%     - The current implementation assumes one oocyte and one cluster.
%     - The repository folder is added to the MATLAB path so the scripts
%       can be run from any working directory.
%
% By Naghmeh Akhavan March 2025

% Ensure all repository scripts are available on the MATLAB path.
script_folder = fileparts(mfilename('fullpath'));
addpath(genpath(script_folder));

% Initialize output structs.
params = struct();
dof = struct();
fd = struct();

%% Cell Geometry And Target Volumes
% Nurse cells.
params.cell_radi_nurse = [0.0005, 0.0005, 0.0005, 0.0005, 0.0005, 0.0005];
params.x_centers_nurse = [1.5, 2.3, 3.1, 1.4, 2.3, 3.0];
params.y_centers_nurse = [1.6, 1.9, 2.0, 3.3, 3.1, 3.0];
params.target_volumes_nurse = [0.5, 0.6, 1.0, 0.6, 0.8, 1.1];

% Oocyte.
params.cell_radi_oocyte = 0.0005;
params.x_centers_oocyte = 3.4;
params.y_centers_oocyte = 2.5;
params.target_volumes_oocyte = 1.95;

% Cluster.
params.cell_radi_cluster = 0.0005;
params.x_centers_cluster = 0.95;
params.y_centers_cluster = 2.5;
params.target_volumes_cluster = 0.2;

% Validate the geometric inputs before building derived arrays.
assert(isequal( ...
    numel(params.cell_radi_nurse), ...
    numel(params.x_centers_nurse), ...
    numel(params.y_centers_nurse), ...
    numel(params.target_volumes_nurse)), ...
    'Nurse-cell parameter vectors must all have the same length.');
assert(isscalar(params.cell_radi_oocyte) && isscalar(params.x_centers_oocyte) && ...
       isscalar(params.y_centers_oocyte) && isscalar(params.target_volumes_oocyte), ...
       'The current solver expects exactly one oocyte.');
assert(isscalar(params.cell_radi_cluster) && isscalar(params.x_centers_cluster) && ...
       isscalar(params.y_centers_cluster) && isscalar(params.target_volumes_cluster), ...
       'The current solver expects exactly one cluster.');

%% Interaction Parameters
% Diffusion coefficients for [nurse cells; oocyte; cluster].
params.dval = [0.001; 0.001; 0.0005];

% Phase-field interaction coefficients.
params.alpha = 100;
params.beta_s = 0.9;     % Epithelial repulsion
params.eta_s = 0.007;    % Epithelial adhesion
params.gamma = 0.01;     % Self-adhesion regularization

% Interaction matrices are ordered as [nurse, oocyte, cluster].
params.beta = [0.25, 0.25, 0.25; ...
               0.25, 0.00, 0.30; ...
               0.25, 0.30, 0.00];
params.eta = [0.003, 0.004, 0.008; ...
              0.004, 0.000, 0.005; ...
              0.008, 0.005, 0.000];

%% Spatial And Temporal Discretization
params.x_length = 5.0;
params.y_length = 5.0;
params.h_space = 0.05;
params.h_time = 0.05;
params.t_final = 500.0;
params.tol = 1e-12;

% Number of saved states for the steady-state and concentration runs.
dof.n_saved = 200;
dof.n_saved_c = 200;

%% Runtime Flags
params.run_steadyQ = 0;   % Run solve_pfm.m
params.run_concenQ = 1;   % Run solve_pfm_concen.m
params.run_tensionQ = 1;  % Label vector-field output as interfacial tension

%% Combined Cell Arrays
% The combined arrays follow the ordering used throughout the solvers:
% nurse cells first, then one oocyte, then one cluster.
params.cell_radi = [params.cell_radi_nurse, params.cell_radi_oocyte, params.cell_radi_cluster];
params.x_centers = [params.x_centers_nurse, params.x_centers_oocyte, params.x_centers_cluster];
params.y_centers = [params.y_centers_nurse, params.y_centers_oocyte, params.y_centers_cluster];
params.target_volumes = [params.target_volumes_nurse, params.target_volumes_oocyte, params.target_volumes_cluster];

dof.n_cells_nurse = numel(params.cell_radi_nurse);
dof.n_cells = numel(params.cell_radi);

%% Concentration-Coupled Parameters
params.h_time_c = 0.05;
params.t_final_c = 1000;

% The concentration field is advanced with smaller internal substeps inside
% each cell update for stability.
params.concen_step = 10;
params.phi0 = 0.5 * ones(1, dof.n_cells);
params.slope2 = 0.02 * ones(1, dof.n_cells);
params.phi0_epi = 0.5;
params.slope2_epi = 0.02;
params.kappa = 0.002;
params.mu_cluster = 0.0;
params.diff_scalar = 0.015;
params.sigma = 0.2;
params.receptor_scalar = 0.005;
params.rhoConcenQ = 1;   % 1: use rho(c) = c, 0: use the nonlinear response

assert(numel(params.phi0) == dof.n_cells, ...
    'params.phi0 must contain one entry per cell.');
assert(numel(params.slope2) == dof.n_cells, ...
    'params.slope2 must contain one entry per cell.');

%% Derived Discretization Metadata
n_space_x = params.x_length / params.h_space + 1;
n_space_y = params.y_length / params.h_space + 1;
steady_steps = params.t_final / params.h_time;
concen_steps = params.t_final_c / params.h_time_c;

assert(abs(n_space_x - round(n_space_x)) < 1e-12 && ...
       abs(n_space_y - round(n_space_y)) < 1e-12, ...
       'x_length and y_length must be integer multiples of h_space.');
assert(abs(n_space_x - n_space_y) < 1e-12, ...
       'The current solver assumes the same number of grid points in x and y.');
assert(abs(steady_steps - round(steady_steps)) < 1e-12, ...
       't_final / h_time must be an integer.');
assert(abs(concen_steps - round(concen_steps)) < 1e-12, ...
       't_final_c / h_time_c must be an integer.');

dof.n_space = round(n_space_x);
dof.n_time = round(steady_steps);

assert(mod(dof.n_time, dof.n_saved) == 0, ...
       'dof.n_saved must divide the number of steady-state time steps.');
assert(mod(round(concen_steps), dof.n_saved_c) == 0, ...
       'dof.n_saved_c must divide the number of concentration time steps.');

%% Mesh, Geometry Mask, And Finite-Difference Operators
fd.x = 0:params.h_space:params.x_length;
fd.y = 0:params.h_space:params.y_length;
[fd.x_mesh, fd.y_mesh] = meshgrid(fd.x, fd.y);

% Fixed epithelial mask used in the adhesion and repulsion terms.
epithelial = tanh(((fd.y_mesh - 2.5).^2 + ((fd.x_mesh - 2.5) * 0.65).^2) ./ 0.5).^100;
fd.epithelial = reshape(epithelial, [dof.n_space^2, 1]);
fd.h_epithelial = 3 * fd.epithelial.^2 - 2 * fd.epithelial.^3;

% Standard five-point Laplacian on the square grid.
id_mat = speye(dof.n_space);
laplacian_1d = (1 / params.h_space^2) * spdiags( ...
    ones(dof.n_space, 1) * [1, 1, -2, 1, 1], ...
    [1 - dof.n_space, -1:1, dof.n_space - 1], ...
    dof.n_space, dof.n_space);
fd.laplacian = kron(laplacian_1d, id_mat) + kron(id_mat, laplacian_1d);
