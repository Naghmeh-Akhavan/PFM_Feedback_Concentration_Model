%GO_SHOW_RESULTS Generate publication-ready plots and videos from saved results.
%   This script loads solver output from the results folder and writes
%   figures or videos to the visuals folder.
%
%   Data sources:
%     results/results.mat
%         Output produced by solve_pfm.m through go_solver.m.
%     results/concentrationsolsaved.mat
%         Output produced by solve_pfm_concen.m through go_solver.m.
%
%   Visualization modes:
%     1  Surface video
%     2  Contour-boundary video
%     3  Filled-contour video
%     4  Residual history plot
%     5  Cluster position and speed plots
%     6  y-slice video for concentration and diffusion coefficient
%     7  Vector-field video
%
%   Field options:
%     1  Cells
%     2  Concentration
%     3  Interfacial tension / receptor force
%     4  Chemical force
%     5  Diffusion coefficient
%
%   Notes:
%     - Modes 4, 5, 6, and 7 do not depend on all field choices.
%     - Modes 2 through 7 require concentration results, except mode 4,
%       which can display steady-state residuals.
%     - Concentration is scaled by params.receptor_scalar by default to
%       preserve the historical visualization convention used in this repo.

%% Visualization Constants
VISUAL_SURFACE = 1;
VISUAL_CONTOUR = 2;
VISUAL_CONTOURF = 3;
VISUAL_RESIDUALS = 4;
VISUAL_CLUSTER_DYNAMICS = 5;
VISUAL_Y_SLICE = 6;
VISUAL_VECTOR_FIELD = 7;

FIELD_CELLS = 1;
FIELD_CONCENTRATION = 2;
FIELD_RECEPTOR = 3;
FIELD_CHEMO = 4;
FIELD_DIFFUSION = 5;

%% User Options
% Choose whether to load steady-state or concentration-coupled results.
use_concentration = 1;

% Select the visualization mode and field type.
type_visual = 2;
type_f = 1;

% Optional display controls.
overlay_vector_fieldQ = 0;           % Only used for contour videos of cells
scale_concentration_by_receptorQ = 1;
slice_y_value = 2.5;
resolution_scale = 1;
frame_rate = 30;

%% Load Saved Results
script_folder = fileparts(mfilename('fullpath'));
results_folder = fullfile(script_folder, 'results');
visuals_folder = fullfile(script_folder, 'visuals');

if ~exist(visuals_folder, 'dir')
    mkdir(visuals_folder);
end

data = load_results_struct(results_folder, use_concentration);
fd = data.fd;
params = data.params;
dof = data.dof;

total_step = get_saved_frame_count(data, use_concentration);
step_time = get_saved_frame_time(params, dof, use_concentration);

if total_step < 1
    error('The selected results file does not contain any saved frames.');
end

%% Dispatch The Requested Visualization
switch type_visual
    case VISUAL_SURFACE
        plot_request = build_plot_request( ...
            data, params, dof, use_concentration, type_f, total_step, scale_concentration_by_receptorQ);
        output_file = fullfile(visuals_folder, sprintf('%s_surface.avi', plot_request.file_label));
        write_surface_video(plot_request, fd, dof, step_time, frame_rate, resolution_scale, output_file);

    case VISUAL_CONTOUR
        plot_request = build_plot_request( ...
            data, params, dof, use_concentration, type_f, total_step, scale_concentration_by_receptorQ);
        output_file = fullfile(visuals_folder, sprintf('%s_contour_bounds.avi', plot_request.file_label));
        write_contour_video( ...
            plot_request, fd, dof, step_time, frame_rate, resolution_scale, output_file, ...
            overlay_vector_fieldQ, data, use_concentration);

    case VISUAL_CONTOURF
        plot_request = build_plot_request( ...
            data, params, dof, use_concentration, type_f, total_step, scale_concentration_by_receptorQ);
        output_file = fullfile(visuals_folder, sprintf('%s_contourf.avi', plot_request.file_label));
        write_filled_contour_video(plot_request, fd, dof, step_time, frame_rate, resolution_scale, output_file);

    case VISUAL_RESIDUALS
        residuals = get_residual_history(data, use_concentration);
        if use_concentration
            output_file = fullfile(visuals_folder, 'concentration_residual_history.png');
            plot_title = 'Concentration-Coupled Residual History';
        else
            output_file = fullfile(visuals_folder, 'steady_state_residual_history.png');
            plot_title = 'Steady-State Residual History';
        end
        plot_residual_history(residuals, step_time, plot_title, output_file);

    case VISUAL_CLUSTER_DYNAMICS
        require_concentration_results(use_concentration, type_visual);
        cluster_history = squeeze(data.sol_cells_saved_c(:, :, dof.n_cells, 1:total_step));
        write_cluster_dynamics_plots(cluster_history, fd.x, step_time, visuals_folder);

    case VISUAL_Y_SLICE
        require_concentration_results(use_concentration, type_visual);
        output_file = fullfile(visuals_folder, 'concentration_diffusion_y_slice.avi');
        write_y_slice_video( ...
            data.concen(:, :, 1:total_step), data.diff_coef(:, :, 1:total_step), ...
            fd, slice_y_value, step_time, frame_rate, resolution_scale, output_file);

    case VISUAL_VECTOR_FIELD
        require_concentration_results(use_concentration, type_visual);
        if params.run_tensionQ
            force_name = 'Interfacial Tension Force';
        else
            force_name = 'Chemical Force';
        end
        output_file = fullfile( ...
            visuals_folder, sprintf('vector_field_%s.avi', lower(strrep(force_name, ' ', '_'))));
        write_vector_field_video( ...
            data.tension_x(:, :, 1:total_step), data.tension_y(:, :, 1:total_step), ...
            fd, force_name, step_time, frame_rate, resolution_scale, output_file);

    otherwise
        error('Invalid type_visual option. Choose 1, 2, 3, 4, 5, 6, or 7.');
end

fprintf('Saved visualization output to %s\n', visuals_folder);

function data = load_results_struct(results_folder, use_concentration)
%LOAD_RESULTS_STRUCT Load the requested results file and validate fields.

    if use_concentration
        results_file = fullfile(results_folder, 'concentrationsolsaved.mat');
        required_fields = {'fd', 'params', 'dof', 'sol_cells_saved_c', 'residuals_c', ...
                           'concen', 'diff_coef', 'receptor', 'chemo', 'tension_x', 'tension_y'};
    else
        results_file = fullfile(results_folder, 'results.mat');
        required_fields = {'fd', 'params', 'dof', 'sol_cells_saved', 'residuals'};
    end

    if ~exist(results_file, 'file')
        error('Results file not found: %s', results_file);
    end

    data = load(results_file);
    missing_fields = required_fields(~isfield(data, required_fields));
    if ~isempty(missing_fields)
        error('The results file %s is missing required variables: %s', ...
            results_file, strjoin(missing_fields, ', '));
    end
end

function total_step = get_saved_frame_count(data, use_concentration)
%GET_SAVED_FRAME_COUNT Determine how many saved frames contain valid data.

    residuals = get_residual_history(data, use_concentration);
    if ~isempty(residuals)
        total_step = numel(residuals);
        return
    end

    if use_concentration
        total_step = size(data.concen, 3);
    else
        total_step = size(data.sol_cells_saved, 4);
    end
end

function step_time = get_saved_frame_time(params, dof, use_concentration)
%GET_SAVED_FRAME_TIME Physical time between saved frames.

    if use_concentration
        step_time = params.t_final_c / dof.n_saved_c;
    else
        step_time = params.t_final / dof.n_saved;
    end
end

function residuals = get_residual_history(data, use_concentration)
%GET_RESIDUAL_HISTORY Return the residual vector associated with the run.

    if use_concentration
        residuals = data.residuals_c(:).';
    else
        residuals = data.residuals(:).';
    end
end

function require_concentration_results(use_concentration, type_visual)
%REQUIRE_CONCENTRATION_RESULTS Guard modes that need concentration output.

    if ~use_concentration
        error('Visualization mode %d requires concentration results. Set use_concentration = 1.', type_visual);
    end
end

function plot_request = build_plot_request(data, params, dof, use_concentration, type_f, total_step, scale_concentration_by_receptorQ)
%BUILD_PLOT_REQUEST Select the requested field and attach metadata.

    plot_request = struct();

    switch type_f
        case 1
            if use_concentration
                plot_request.field = data.sol_cells_saved_c(:, :, :, 1:total_step);
            else
                plot_request.field = data.sol_cells_saved(:, :, :, 1:total_step);
            end
            plot_request.display_name = 'Cells';
            plot_request.file_label = 'cells';
            plot_request.is_cells = true;
            plot_request.clim = [];

        case 2
            require_concentration_results(use_concentration, type_f);
            if scale_concentration_by_receptorQ
                plot_request.field = params.receptor_scalar * data.concen(:, :, 1:total_step);
                plot_request.display_name = 'Scaled Concentration';
            else
                plot_request.field = data.concen(:, :, 1:total_step);
                plot_request.display_name = 'Concentration';
            end
            plot_request.file_label = 'concentration';
            plot_request.is_cells = false;
            plot_request.clim = expanded_limits(plot_request.field);

        case 3
            require_concentration_results(use_concentration, type_f);
            plot_request.field = data.receptor(:, :, 1:total_step);
            plot_request.display_name = 'Interfacial Tension';
            plot_request.file_label = 'interfacial_tension';
            plot_request.is_cells = false;
            plot_request.clim = expanded_limits(plot_request.field);

        case 4
            require_concentration_results(use_concentration, type_f);
            plot_request.field = data.chemo(:, :, 1:total_step);
            plot_request.display_name = 'Chemical Force';
            plot_request.file_label = 'chemical_force';
            plot_request.is_cells = false;
            plot_request.clim = expanded_limits(plot_request.field);

        case 5
            require_concentration_results(use_concentration, type_f);
            plot_request.field = data.diff_coef(:, :, 1:total_step);
            plot_request.display_name = 'Diffusion Coefficient';
            plot_request.file_label = 'diffusion_coefficient';
            plot_request.is_cells = false;
            plot_request.clim = expanded_limits(plot_request.field);

        otherwise
            error('Invalid type_f option. Choose 1, 2, 3, 4, or 5.');
    end

    plot_request.nurse_count = dof.n_cells_nurse;
end

function write_surface_video(plot_request, fd, dof, step_time, frame_rate, resolution_scale, output_file)
%WRITE_SURFACE_VIDEO Render a 3-D surface video.

    writer = VideoWriter(output_file);
    writer.FrameRate = frame_rate;
    open(writer);

    fig = figure('Units', 'pixels', 'Position', [0, 0, 560 * resolution_scale, 420 * resolution_scale]);
    frame_total = get_frame_total(plot_request.field);

    for frame_idx = 1:frame_total
        clf(fig);
        ax = axes(fig);
        hold(ax, 'on');

        if plot_request.is_cells
            draw_cell_surface_frame(ax, fd, dof, plot_request.field(:, :, :, frame_idx), plot_request.nurse_count);
        else
            surf(ax, fd.x, fd.y, plot_request.field(:, :, frame_idx), 'EdgeColor', 'none');
            view(ax, -50, 60);
            xlim(ax, [fd.x(1), fd.x(end)]);
            ylim(ax, [fd.y(1), fd.y(end)]);
            axis(ax, 'tight');
            colorbar(ax);
            clim(ax, plot_request.clim);
        end

        xlabel(ax, 'x');
        ylabel(ax, 'y');
        title(ax, sprintf('%s at Time: %.2f', plot_request.display_name, (frame_idx - 1) * step_time));
        drawnow;
        writeVideo(writer, getframe(fig));
    end

    close(writer);
    close(fig);
end

function write_contour_video(plot_request, fd, dof, step_time, frame_rate, resolution_scale, output_file, overlay_vector_fieldQ, data, use_concentration)
%WRITE_CONTOUR_VIDEO Render a contour-boundary video.

    writer = VideoWriter(output_file);
    writer.FrameRate = frame_rate;
    open(writer);

    fig = figure('Units', 'pixels', 'Position', [0, 0, 560 * resolution_scale, 420 * resolution_scale]);
    frame_total = get_frame_total(plot_request.field);

    for frame_idx = 1:frame_total
        clf(fig);
        ax = axes(fig);
        hold(ax, 'on');

        if plot_request.is_cells
            draw_cell_contour_frame(ax, fd, dof, plot_request.field(:, :, :, frame_idx), plot_request.nurse_count);

            if overlay_vector_fieldQ && use_concentration
                quiver(ax, fd.x, fd.y, data.tension_x(:, :, frame_idx), data.tension_y(:, :, frame_idx), ...
                    'Color', [1, 0, 1], 'LineWidth', 1.0);
            end
        else
            contour(ax, fd.x, fd.y, plot_request.field(:, :, frame_idx), 20, 'LineWidth', 1.2);
            colorbar(ax);
            clim(ax, plot_request.clim);
        end

        xlabel(ax, 'x');
        ylabel(ax, 'y');
        axis(ax, 'equal');
        axis(ax, 'tight');
        title(ax, sprintf('%s at Time: %.2f', plot_request.display_name, (frame_idx - 1) * step_time));
        drawnow;
        writeVideo(writer, getframe(fig));
    end

    close(writer);
    close(fig);
end

function write_filled_contour_video(plot_request, fd, dof, step_time, frame_rate, resolution_scale, output_file)
%WRITE_FILLED_CONTOUR_VIDEO Render a filled-contour video.

    writer = VideoWriter(output_file);
    writer.FrameRate = frame_rate;
    open(writer);

    fig = figure('Units', 'pixels', 'Position', [0, 0, 560 * resolution_scale, 420 * resolution_scale]);
    frame_total = get_frame_total(plot_request.field);

    for frame_idx = 1:frame_total
        clf(fig);
        ax = axes(fig);
        hold(ax, 'on');

        if plot_request.is_cells
            draw_cell_filled_contour_frame(ax, fd, dof, plot_request.field(:, :, :, frame_idx), plot_request.nurse_count);
            colormap(ax, flipud(gray));
        else
            contourf(ax, fd.x, fd.y, plot_request.field(:, :, frame_idx), 20, 'LineStyle', 'none');
            colorbar(ax);
            clim(ax, plot_request.clim);
        end

        xlabel(ax, 'x');
        ylabel(ax, 'y');
        axis(ax, 'equal');
        axis(ax, 'tight');
        title(ax, sprintf('%s at Time: %.2f', plot_request.display_name, (frame_idx - 1) * step_time));
        drawnow;
        writeVideo(writer, getframe(fig));
    end

    close(writer);
    close(fig);
end

function plot_residual_history(residuals, step_time, plot_title, output_file)
%PLOT_RESIDUAL_HISTORY Save a residual-history figure.

    time = (1:numel(residuals)) * step_time;
    fig = figure('Units', 'pixels', 'Position', [0, 0, 700, 500]);
    plot(time, residuals, 'LineWidth', 2);
    xlabel('Time');
    ylabel('Relative Residual');
    title(plot_title);
    grid on;
    saveas(fig, output_file);
    close(fig);
end

function write_cluster_dynamics_plots(cluster_history, x_coordinates, step_time, visuals_folder)
%WRITE_CLUSTER_DYNAMICS_PLOTS Save cluster centroid and speed plots.

    frame_total = size(cluster_history, 3);
    avg_x_positions = nan(1, frame_total);

    for frame_idx = 1:frame_total
        column_weights = sum(cluster_history(:, :, frame_idx), 1);
        total_weight = sum(column_weights);
        if total_weight > 0
            avg_x_positions(frame_idx) = sum(column_weights .* x_coordinates) / total_weight;
        end
    end

    speeds = diff(avg_x_positions) / step_time;
    time = (0:frame_total - 1) * step_time;

    fig = figure('Units', 'pixels', 'Position', [0, 0, 700, 500]);
    plot(time, avg_x_positions, 'LineWidth', 3);
    xlabel('Time');
    ylabel('Average X Position');
    title('Weighted Average X Position Of The Cluster Over Time');
    grid on;
    saveas(fig, fullfile(visuals_folder, 'cluster_avg_x_position.png'));
    close(fig);

    fig = figure('Units', 'pixels', 'Position', [0, 0, 700, 500]);
    plot(time(2:end), speeds, 'LineWidth', 3);
    xlabel('Time');
    ylabel('Speed');
    title('Cluster Speed Over Time');
    grid on;
    saveas(fig, fullfile(visuals_folder, 'cluster_speed.png'));
    close(fig);

    fig = figure('Units', 'pixels', 'Position', [0, 0, 700, 500]);
    yyaxis left
    plot(time, avg_x_positions, 'b-', 'LineWidth', 3);
    ylabel('Average X Position');

    yyaxis right
    plot(time(2:end), speeds, 'r-', 'LineWidth', 3);
    ylabel('Speed');

    xlabel('Time');
    title('Cluster Position And Speed Over Time');
    grid on;
    legend('Average X Position', 'Speed', 'Location', 'best');
    saveas(fig, fullfile(visuals_folder, 'cluster_position_speed.png'));
    close(fig);
end

function write_y_slice_video(concen, diff_coef, fd, y_value, step_time, frame_rate, resolution_scale, output_file)
%WRITE_Y_SLICE_VIDEO Animate a horizontal slice of concentration and diffusion.

    [~, y_slice_index] = min(abs(fd.y - y_value));
    y_slice_value = fd.y(y_slice_index);

    conc_limits = expanded_limits(concen(y_slice_index, :, :));
    diff_limits = expanded_limits(diff_coef(y_slice_index, :, :));
    frame_total = size(concen, 3);

    writer = VideoWriter(output_file);
    writer.FrameRate = frame_rate;
    open(writer);

    fig = figure('Units', 'pixels', 'Position', [0, 0, 560 * resolution_scale, 420 * resolution_scale]);
    for frame_idx = 1:frame_total
        clf(fig);

        yyaxis left
        plot(fd.x, squeeze(concen(y_slice_index, :, frame_idx)), 'b-', 'LineWidth', 3);
        ylabel('Concentration');
        ylim(conc_limits);
        set(gca, 'YColor', 'b');

        yyaxis right
        plot(fd.x, squeeze(diff_coef(y_slice_index, :, frame_idx)), 'g-', 'LineWidth', 3);
        ylabel('Diffusion Coefficient');
        ylim(diff_limits);
        set(gca, 'YColor', 'g');

        xlabel('x');
        title(sprintf('Concentration And Diffusion At y = %.2f, Time = %.2f', ...
            y_slice_value, (frame_idx - 1) * step_time));
        legend('Concentration', 'Diffusion Coefficient', 'Location', 'best');
        grid on;
        drawnow;
        writeVideo(writer, getframe(fig));
    end

    close(writer);
    close(fig);
end

function write_vector_field_video(tension_x, tension_y, fd, force_name, step_time, frame_rate, resolution_scale, output_file)
%WRITE_VECTOR_FIELD_VIDEO Animate the stored vector field.

    frame_total = size(tension_x, 3);
    writer = VideoWriter(output_file);
    writer.FrameRate = frame_rate;
    open(writer);

    fig = figure('Units', 'pixels', 'Position', [0, 0, 560 * resolution_scale, 420 * resolution_scale]);
    for frame_idx = 1:frame_total
        clf(fig);
        ax = axes(fig);
        quiver(ax, fd.x, fd.y, tension_x(:, :, frame_idx), tension_y(:, :, frame_idx), ...
            'AutoScale', 'on', 'AutoScaleFactor', 2, 'LineWidth', 1.0);
        xlabel(ax, 'x');
        ylabel(ax, 'y');
        axis(ax, 'equal');
        axis(ax, 'tight');
        title(ax, sprintf('%s at Time: %.2f', force_name, (frame_idx - 1) * step_time));
        drawnow;
        writeVideo(writer, getframe(fig));
    end

    close(writer);
    close(fig);
end

function draw_cell_surface_frame(ax, fd, dof, cell_frame, nurse_count)
%DRAW_CELL_SURFACE_FRAME Draw all cells as semi-transparent surfaces.

    epithelial_surface = reshape(fd.epithelial, dof.n_space, dof.n_space);
    surf(ax, fd.x, fd.y, epithelial_surface, ...
        'FaceColor', [0.25, 0.25, 0.25], 'FaceAlpha', 0.30, 'EdgeColor', 'none');

    cell_count = size(cell_frame, 3);
    for cell_idx = 1:cell_count
        surf(ax, fd.x, fd.y, cell_frame(:, :, cell_idx), ...
            'FaceColor', get_cell_color(cell_idx, cell_count, nurse_count), ...
            'FaceAlpha', 0.65, 'EdgeColor', 'none');
    end

    shading(ax, 'interp');
    xlim(ax, [fd.x(1), fd.x(end)]);
    ylim(ax, [fd.y(1), fd.y(end)]);
    zlim(ax, [0, 1]);
    view(ax, -50, 60);
end

function draw_cell_contour_frame(ax, fd, dof, cell_frame, nurse_count)
%DRAW_CELL_CONTOUR_FRAME Draw cell boundaries together with the epithelial mask.

    contour(ax, fd.x, fd.y, reshape(fd.epithelial, dof.n_space, dof.n_space), ...
        [0.5, 0.5], 'k', 'LineWidth', 2.5);

    cell_count = size(cell_frame, 3);
    for cell_idx = 1:cell_count
        contour(ax, fd.x, fd.y, cell_frame(:, :, cell_idx), [0.5, 0.5], ...
            'LineColor', get_cell_color(cell_idx, cell_count, nurse_count), 'LineWidth', 1.5);
    end
end

function draw_cell_filled_contour_frame(ax, fd, dof, cell_frame, nurse_count)
%DRAW_CELL_FILLED_CONTOUR_FRAME Draw filled cell contours and the epithelial mask.

    cell_count = size(cell_frame, 3);
    for cell_idx = 1:cell_count
        contourf(ax, fd.x, fd.y, cell_frame(:, :, cell_idx), 0.5:0.1:1.0, 'LineStyle', 'none');
    end

    contour(ax, fd.x, fd.y, reshape(fd.epithelial, dof.n_space, dof.n_space), ...
        [0.5, 0.5], 'k', 'LineWidth', 2.5);

    for cell_idx = 1:cell_count
        contour(ax, fd.x, fd.y, cell_frame(:, :, cell_idx), [0.5, 0.5], ...
            'LineColor', get_cell_color(cell_idx, cell_count, nurse_count), 'LineWidth', 1.0);
    end
end

function color = get_cell_color(cell_idx, cell_count, nurse_count)
%GET_CELL_COLOR Assign consistent colors by cell type.

    if cell_idx <= nurse_count
        color = [0.85, 0.15, 0.15];
    elseif cell_idx == cell_count - 1
        color = [0.15, 0.35, 0.85];
    else
        color = [0.10, 0.65, 0.20];
    end
end

function limits = expanded_limits(values)
%EXPANDED_LIMITS Return stable plotting limits for possibly constant data.

    finite_values = values(isfinite(values));
    if isempty(finite_values)
        limits = [0, 1];
        return
    end

    data_min = min(finite_values);
    data_max = max(finite_values);

    if data_min == data_max
        pad = max(1, abs(data_min)) * 0.1;
        limits = [data_min - pad, data_max + pad];
    else
        limits = [data_min, data_max];
    end
end

function frame_total = get_frame_total(field_data)
%GET_FRAME_TOTAL Return the number of saved frames in a field array.

    if ndims(field_data) == 4
        frame_total = size(field_data, 4);
    else
        frame_total = size(field_data, 3);
    end
end
