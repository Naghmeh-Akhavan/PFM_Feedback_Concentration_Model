%GO_SOLVER Run the steady-state and concentration solvers for the PFM model.
%   This script is the main entry point
%   1. solve the phase-field model without concentration to obtain a
%      steady-state cell configuration, and
%   2. restart from that saved state to solve the concentration-coupled
%      model.
%
%   Runtime control is handled through the boolean flags in config.m:
%     params.run_steadyQ  Run the steady-state solver.
%     params.run_concenQ  Run the concentration solver.
%
%   Output files:
%     results/results.mat
%         Saved output from solve_pfm.m.
%     results/concentrationsolsaved.mat
%         Saved output from solve_pfm_concen.m.
%
%   Notes:
%     - The concentration solve expects a previously saved steady-state
%       solution.
%     - The penultimate saved steady-state frame is used as the initial
%       condition for the concentration run, matching the original script
%       workflow.

clear variables

% Resolve output locations relative to this script so the workflow is
% independent of the current working directory.
script_folder = fileparts(mfilename('fullpath'));
results_folder = fullfile(script_folder, 'results');
steady_results_file = fullfile(results_folder, 'results.mat');
concen_results_file = fullfile(results_folder, 'concentrationsolsaved.mat');

% Create the output folder if it does not already exist.
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

% Load the model configuration and optionally run the steady-state solve.
[fd, params, dof] = config();
if params.run_steadyQ
    tstart = tic;
    [sol_cells_saved, residuals] = solve_pfm(fd, params, dof);
    save(steady_results_file);
    fprintf('Time to solve the steady-state problem: %.3f seconds.\n', toc(tstart));
end

% Load the steady-state solution needed to initialize the concentration
% solver.
if exist(steady_results_file, 'file')
    load(steady_results_file);
else
    error(['No steady-state results file was found at %s.\n' ...
           'Run the steady-state solver first by setting params.run_steadyQ = 1 in config.m.'], ...
          steady_results_file);
end

% Reload the configuration so the concentration-specific time-stepping
% parameters are freshly available before the restart solve.
[fd, params, dof] = config();
if params.run_concenQ
    tstart = tic;

    % Use the penultimate saved steady-state frame as the restart state.
    cell_state = sol_cells_saved(:, :, :, end - 1);
    [sol_cells_saved_c, residuals_c, concen, diff_coef, receptor, chemo, tension_x, tension_y] = ...
        solve_pfm_concen(fd, params, dof, cell_state);

    save(concen_results_file);
    fprintf('Time to solve the concentration-coupled problem: %.3f seconds.\n', toc(tstart));
end
