function u_initial = initialize_u(x_mesh, y_mesh,x_center, y_center, d_zero)
% Generates Gaussian function over given meshgrid
% Inputs:
%   x_mesh, y_mesh          meshgrid for the domain
%   x_center, y_center      center of the Gaussian function
%   d_zero                  standard deviation of the Gaussian function
% Outputs:
%   u_initial               Gaussian function over the meshgrid
%
% By Naghmeh Akhavan, March 2025


tan_input = ((x_mesh-x_center).^2 + (y_mesh-y_center).^2)./(2*sqrt(2*d_zero));
u_initial = 0.5*(1 - tanh(tan_input));

end