% Pseudo-boundary-condition diffusion demo script
% - solves a 1D diffusion–decay PDE with moving boundary-like phi, via pdepe
% - includes variable effective diffusion D(phi), source and secretion terms
% - outputs time-resolved c(x,t), phi(x,t), and D_phi(x,t) with video
% - parameters (D0, phi0, s1, s2, k, v, x0, sigma, R) are provided in params

% Define spatial domain
L = 10;  % Length of the domain
x = linspace(0, L, 200); % Spatial grid
t = linspace(0,200, 400); % Time grid
% Create a meshgrid from x and t
[X, T] = meshgrid(x, t);

% Model parameters
params.D0 = 1;                % Base diffusion coefficient
params.phi0 = 0.0;            % Phase-field threshold
params.s1 = 0.2;              % Gaussian width (unused)
params.s2 = 0.02;             % Sigmoid steepness for D(phi)
params.k = 0.1;               % Decay rate
params.v = 0.005;             % Boundary velocity
params.x0 = 5;                % Initial boundary position
params.sigma = 0.2;           % Secretion rate
params.R = 1;                 % Boundary half-width
pdefun = @(x, t, c, c_x) pde_function(x, t, c, c_x, params);
icfun = @(x) initial_condition(x, params);
bcfun = @(xl, cl, xr, cr, t) boundary_condition(xl, cl, xr, cr, t, params);

% Solve PDE using MATLAB's pdepe solver (method of lines)
sol = pdepe(0, pdefun, icfun, bcfun, x, t);
c = sol(:,:,1);  % Extract concentration field

% Save solution for later analysis
save('pde_concen_solution.mat');

% Reconstruct moving phase field and diffusion coefficient
[X, T] = meshgrid(x, t);
x0 = params.x0 + params.v * T;  % Time-dependent boundary position
R = params.R;
phi = 0.5 * (tanh((X - (x0 - R)) / 0.1) - tanh((X - (x0 + R)) / 0.1));

% Generate video output showing solution evolution
% Extract and cache parameters for video loop
D0 = params.D0;
phi0 = params.phi0;
s1 = params.s1;
s2 = params.s2;
k = params.k;
v = params.v;
x0 = params.x0 + v*T;
R = params.R;

% Compute effective diffusion as sigmoidal function of phase field
phi = 0.5 * (tanh((x - (x0 - R)) / 0.1) - tanh((x - (x0 + R)) / 0.1));
D_phi = D0 ./ ((1 + exp((phi - phi0)./s2))); 
v = VideoWriter('pde_solution.avi');
open(v);
p = plot(x, c(1,:), 'b', x, phi(1,:), 'r--', x, D_phi(1,:), 'g-.',"LineWidth",2);
ylim([0 max(c(end,:))]);
title(sprintf('Time: %.2f', t(1)));
legend('c(x,t)', 'phi(x,t)', 'D\_phi(x,t)');
frame = getframe(gcf);
writeVideo(v, frame);
% Animate solution over all time steps
for i = 1:length(t)
    % Update plot data for concentration, phase field, and diffusion
    set(p(1), 'YData', c(i,:));      % c(x,t)
    set(p(2), 'YData', phi(i,:));    % phi(x,t)
    set(p(3), 'YData', D_phi(i,:));  % D(phi(x,t))
    title(sprintf('Time: %.2f', t(i)));
    frame = getframe(gcf);
    writeVideo(v, frame);
    pause(0.1);
end

% Finalize video file
close(v);

% Display snapshot at first time step for reference
t_step = 1; 
plot(x, c(t_step,:), 'b', x, phi(t_step,:), 'r--', x, D_phi(t_step,:), 'g-.',"LineWidth",2);
ylim([0 max(c(end,:))]);
title(sprintf('Time: %.2f', t(t_step)));
legend('c(x,t)', 'phi(x,t)', 'D\_phi(x,t)');

% PDE residual function for pdepe solver
% Solves: c_t = d/dx[D(phi) * dc/dx] - k*c
function [c_t, f, s] = pde_function(x, t, c, c_x, params)
    % Extract parameters
    D0 = params.D0;
    phi0 = params.phi0;
    s2 = params.s2;
    k = params.k;
    v = params.v;
    x0 = params.x0 + v*t;
    R = params.R;

    % Compute phase field (sharp boundary profile)
    phi = 0.5 * (tanh((x - (x0 - R)) / 0.1) - tanh((x - (x0 + R)) / 0.1));

    % Compute effective diffusion coefficient (sigmoid in phase field)
    D_phi = D0 / ((1 + exp((phi - phi0)/s2)));

    % PDE coefficients for pdepe format: c_t = d/dx[f] + s
    c_t = 1;           % Time derivative coefficient
    f = D_phi * c_x;   % Divergence of flux D(phi)*grad(c)
    s = -k * c;        % Decay source term
end

% Initial condition: homogeneous zero concentration
function c0 = initial_condition(x, params)
    c0 = 0 * x;  % Zero concentration throughout domain at t=0
end

% Boundary conditions: Dirichlet form [p + q*c_x = 0]
% Left (x=0): no-flux; Right (x=L): secretion
function [pl, ql, pr, qr] = boundary_condition(xl, cl, xr, cr, t, params)
    % Extract parameters
    D0 = params.D0;
    phi0 = params.phi0;
    s2 = params.s2;
    sigma = params.sigma;  % Secretion rate
    x0 = params.x0 + params.v*t;
    R = params.R;

    % Compute phase field and effective diffusion at left boundary (x=0)
    phiL = 0.5*(tanh(0-(x0-R)/0.1) - tanh(0-(x0+R)/0.1));
    DL = D0 / ((1 + exp((phiL - phi0)/s2)));

    % Compute phase field and effective diffusion at right boundary (x=L=10)
    phiR = 0.5*(tanh(10-(x0-R)/0.1) - tanh(10-(x0+R)/0.1));
    DR = D0 / ((1 + exp((phiR - phi0)/s2)));

    % Left boundary: no-flux condition (p=0, q=1 means c_x=0)
    pl = 0;
    ql = 1;

    % Right boundary: source/secretion at rate sigma
    pr = -sigma / DR;
    qr = 1;
end
