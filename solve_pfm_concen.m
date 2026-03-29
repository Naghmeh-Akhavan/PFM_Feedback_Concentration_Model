function [sol_cells_saved, residuals, concen_saved, diff_coef_saved, receptor_saved, chemo_saved, tension_x_saved, tension_y_saved, concen_iters_saved] = solve_pfm_concen(fd, params, dof, cell_state)
%SOLVE_PFM_CONCEN Advance the phase-field model with concentration dynamics.
%   [sol_cells_saved, residuals, concen_saved, diff_coef_saved, ...
%    receptor_saved, chemo_saved, tension_x_saved, tension_y_saved, ...
%    concen_iters_saved] = solve_pfm_concen(fd, params, dof, cell_state)
%   starts from a steady-state cell configuration and advances the coupled
%   cell, concentration, receptor, and chemotaxis fields.
%
%   Inputs:
%     fd         Finite-difference operators and geometry data.
%     params     Model parameters, including concentration time-stepping
%                settings.
%     dof        Discretization sizes and output cadence.
%     cell_state Steady-state snapshot of size
%                n_space-by-n_space-by-n_cells.
%
%   Outputs:
%     sol_cells_saved    Saved cell phase fields.
%     residuals          Relative cluster residual at saved times.
%     concen_saved       Saved concentration field.
%     diff_coef_saved    Saved effective diffusion coefficient field.
%     receptor_saved     Saved receptor-force field.
%     chemo_saved        Saved chemotactic forcing field.
%     tension_x_saved    Saved x-component of the tension field.
%     tension_y_saved    Saved y-component of the tension field.
%     concen_iters_saved Number of concentration substeps per saved frame.
%
% By Naghmeh Akhavan March 2025

% Override the generic time-stepping values with the concentration run.
dof.n_saved = dof.n_saved_c;
params.h_time = params.h_time_c;
params.t_final = params.t_final_c;
dof.n_time = params.t_final / params.h_time;
save_interval = dof.n_time / dof.n_saved;

% Common indices used throughout the update.
nurse_idx = 1:dof.n_cells_nurse;
oocyte_idx = dof.n_cells - 1;
cluster_idx = dof.n_cells;

% Vectorize the steady-state snapshot for the time-marching loop.
sol_cells = reshape(cell_state, dof.n_space^2, dof.n_cells);

% Preallocate saved outputs.
sol_cells_saved = zeros(dof.n_space, dof.n_space, dof.n_cells, dof.n_saved);
concen_saved = zeros(dof.n_space, dof.n_space, dof.n_saved);
diff_coef_saved = zeros(dof.n_space, dof.n_space, dof.n_saved);
receptor_saved = zeros(dof.n_space, dof.n_space, dof.n_saved);
chemo_saved = zeros(dof.n_space, dof.n_space, dof.n_saved);
tension_x_saved = zeros(dof.n_space, dof.n_space, dof.n_saved);
tension_y_saved = zeros(dof.n_space, dof.n_space, dof.n_saved);
concen_iters_saved = zeros(dof.n_saved, 1);
residuals = zeros(1, dof.n_saved);

% Build first-derivative operators with Neumann boundary conditions.
id_mat = speye(dof.n_space);
central_diff_1d = -(1 / (2 * params.h_space)) * spdiags( ...
    ones(dof.n_space, 1) * [1, -1], [-1, 1], dof.n_space, dof.n_space);
central_diff_1d(1, 1) = -1 / (2 * params.h_space);
central_diff_1d(dof.n_space, dof.n_space) = 1 / (2 * params.h_space);
central_diff_x = kron(central_diff_1d, id_mat);
central_diff_y = kron(id_mat, central_diff_1d);

% Initialize the concentration field.
concen = zeros(dof.n_space^2, 1);

% Cache frequently used parameters to keep the update expressions readable.
phi0 = params.phi0;
slope2 = params.slope2;
phi0_epi = params.phi0_epi;
slope2_epi = params.slope2_epi;
kappa = params.kappa;
mu_cluster = params.mu_cluster;

save_count = 0;
progress_save_interval = 50;
progress_clock = tic;

for i = 1:dof.n_time
    % Suppress tiny round-off drift before building nonlinear terms.
    sol_cells = round(sol_cells, 16);
    sol_cells_old = sol_cells;

    % Smoothed indicator h(u) = 3u^2 - 2u^3 and the aggregate nurse phase.
    h_sol_cells = 3 * sol_cells.^2 - 2 * sol_cells.^3;
    phi_sum = sum(h_sol_cells(:, nurse_idx), 2);
    vol_cells = params.h_space^2 * sum(h_sol_cells, 1);

    % Update nurse cells.
    for m = nurse_idx
        sol_cells(:, m) = sol_cells(:, m) + params.h_time * ( ...
            params.dval(1) * fd.laplacian * sol_cells(:, m) ...
            + sol_cells(:, m) .* (1 - sol_cells(:, m)) .* ( ...
                sol_cells(:, m) - 0.5 + (params.target_volumes(m) - vol_cells(m)) ...
                - params.beta(1, 1) * (phi_sum - h_sol_cells(:, m)) ...
                - params.beta(1, 2) * h_sol_cells(:, oocyte_idx) ...
                - params.beta(1, 3) * h_sol_cells(:, cluster_idx) ...
                - params.beta_s * fd.h_epithelial ...
                + fd.laplacian * ( ...
                    params.eta_s * fd.h_epithelial ...
                    + params.eta(1, 1) * (phi_sum - h_sol_cells(:, m)) ...
                    + params.eta(1, 2) * h_sol_cells(:, oocyte_idx) ...
                    + params.eta(1, 3) * h_sol_cells(:, cluster_idx) ...
                    + params.gamma * h_sol_cells(:, m))));
    end

    % Update the oocyte.
    sol_cells(:, oocyte_idx) = sol_cells(:, oocyte_idx) + params.h_time * ( ...
        params.dval(2) * fd.laplacian * sol_cells(:, oocyte_idx) ...
        + sol_cells(:, oocyte_idx) .* (1 - sol_cells(:, oocyte_idx)) .* ( ...
            sol_cells(:, oocyte_idx) - 0.5 + (params.target_volumes(oocyte_idx) - vol_cells(oocyte_idx)) ...
            - params.beta(2, 1) * phi_sum ...
            - params.beta(2, 3) * h_sol_cells(:, cluster_idx) ...
            - params.beta_s * fd.h_epithelial ...
            + fd.laplacian * ( ...
                params.eta_s * fd.h_epithelial ...
                + params.eta(2, 1) * phi_sum ...
                + params.eta(2, 3) * h_sol_cells(:, cluster_idx) ...
                + params.gamma * h_sol_cells(:, oocyte_idx))));

    % Build the diffusion coefficient from the current phase fields and
    % rescale it so diffusion vanishes at the peak of the last nurse cell.
    exp_term = exp(reshape(fd.h_epithelial, dof.n_space, dof.n_space) - phi0_epi) ./ slope2_epi;
    for k = 1:dof.n_cells
        exp_term = exp_term + exp(reshape(h_sol_cells(:, k), dof.n_space, dof.n_space) - phi0(k)) ./ slope2(k);
    end
    diff_coef = 1 ./ (1 + exp_term);
    diff_coef = reshape(diff_coef, [dof.n_space^2, 1]);

    [~, ref_idx] = max(sol_cells(:, nurse_idx(end)));
    diff_coef = diff_coef - diff_coef(ref_idx);

    max_diff_coef = max(diff_coef);
    if max_diff_coef > 0
        diff_coef = diff_coef / max_diff_coef;
    else
        diff_coef = zeros(size(diff_coef));
    end

    diff_coef = 3 * diff_coef.^2 - 2 * diff_coef.^3;
    diff_coef = params.diff_scalar * diff_coef;

    % Secretion is proportional to the oocyte interface gradient.
    [oocyte_grad_x, oocyte_grad_y] = gradient( ...
        reshape(sol_cells(:, oocyte_idx), dof.n_space, dof.n_space), params.h_space);
    secretion = params.sigma * hypot(oocyte_grad_x, oocyte_grad_y);
    secretion = reshape(secretion, [dof.n_space^2, 1]);

    % Take concentration substeps for stability within each cell update.
    for substep = 1:params.concen_step
        [concen_x, concen_y] = gradient(reshape(concen, [dof.n_space, dof.n_space]), params.h_space);
        concen_x = reshape(concen_x, [dof.n_space^2, 1]);
        concen_y = reshape(concen_y, [dof.n_space^2, 1]);

        delta_concen = (params.h_time / params.concen_step) * ( ...
            central_diff_x * (diff_coef .* concen_x) ...
            + central_diff_y * (diff_coef .* concen_y) ...
            - kappa * concen + secretion);
        concen = concen + delta_concen;
    end

    % Build receptor-driven tension from the concentration and cluster
    % interface geometry.
    [cluster_x, cluster_y] = gradient( ...
        reshape(sol_cells(:, cluster_idx), dof.n_space, dof.n_space), params.h_space);
    cluster_x = reshape(cluster_x, [dof.n_space^2, 1]);
    cluster_y = reshape(cluster_y, [dof.n_space^2, 1]);

    rho_concen = concen.^3 ./ ((concen.^2 + 1) .* (concen + 1));

    nurse_presence = sum(sol_cells(:, nurse_idx), 2);
    turning_sign = sign(-concen_y .* cluster_x + concen_x .* cluster_y);
    tension_x = rho_concen .* sol_cells(:, cluster_idx) .* nurse_presence .* (turning_sign .* cluster_y);
    tension_y = -rho_concen .* sol_cells(:, cluster_idx) .* nurse_presence .* (turning_sign .* cluster_x);
    receptor_force = -params.receptor_scalar * ( ...
        central_diff_x * tension_x + central_diff_y * tension_y);

    % Combine receptor forcing, chemotaxis, and phase-field dynamics for
    % the cluster update.
    chem_x = sol_cells(:, cluster_idx) .* concen_x;
    chem_y = sol_cells(:, cluster_idx) .* concen_y;
    chemotaxis_force = -mu_cluster * (central_diff_x * chem_x + central_diff_y * chem_y);

    sol_cells(:, cluster_idx) = sol_cells(:, cluster_idx) + params.h_time * ( ...
        receptor_force ...
        + chemotaxis_force ...
        + params.dval(3) * fd.laplacian * sol_cells(:, cluster_idx) ...
        + sol_cells(:, cluster_idx) .* (1 - sol_cells(:, cluster_idx)) .* ( ...
            sol_cells(:, cluster_idx) - 0.5 + params.alpha * (params.target_volumes(cluster_idx) - vol_cells(cluster_idx)) ...
            - params.beta(3, 1) * phi_sum ...
            - params.beta(3, 2) * h_sol_cells(:, oocyte_idx) ...
            - params.beta_s * fd.h_epithelial ...
            + fd.laplacian * ( ...
                params.eta_s * fd.h_epithelial ...
                + params.eta(3, 1) * phi_sum ...
                + params.eta(3, 2) * h_sol_cells(:, oocyte_idx) ...
                + params.gamma * h_sol_cells(:, cluster_idx))));

    res = norm(sol_cells(:, cluster_idx) - sol_cells_old(:, cluster_idx)) / ...
        norm(sol_cells_old(:, cluster_idx));

    if mod(i, save_interval) == 0
        save_count = save_count + 1;

        sol_cells_saved(:, :, :, save_count) = reshape(sol_cells, dof.n_space, dof.n_space, dof.n_cells);
        concen_saved(:, :, save_count) = reshape(concen, dof.n_space, dof.n_space);
        diff_coef_saved(:, :, save_count) = reshape(diff_coef, dof.n_space, dof.n_space);
        receptor_saved(:, :, save_count) = reshape(receptor_force, dof.n_space, dof.n_space);
        tension_x_saved(:, :, save_count) = reshape(tension_x, dof.n_space, dof.n_space);
        tension_y_saved(:, :, save_count) = reshape(tension_y, dof.n_space, dof.n_space);
        chemo_saved(:, :, save_count) = reshape(chemotaxis_force, dof.n_space, dof.n_space);
        concen_iters_saved(save_count) = params.concen_step;
        residuals(save_count) = res;

        fprintf(['Saved time step %d of %d. Relative cluster residual: %.3e. ' ...
                 'Previous cluster tip value: %.3e.\n'], ...
                i, dof.n_time, res, max(sol_cells_old(:, cluster_idx)));

        if mod(save_count, progress_save_interval) == 0
            elapsed = toc(progress_clock);
            completed_percent = 100 * save_count / dof.n_saved;
            remaining_blocks = (dof.n_saved - save_count) / progress_save_interval;
            fprintf('Completed %.0f%% of saved outputs. Estimated time remaining: %.1f seconds.\n', ...
                completed_percent, elapsed * remaining_blocks);
            progress_clock = tic;
        end

        if res < params.tol
            break
        end
    end
end

residuals = residuals(1:save_count);
concen_iters_saved = concen_iters_saved(1:save_count);
