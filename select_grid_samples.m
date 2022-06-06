function x_grid_sample = select_grid_samples(n_samples_per_dim,x_lower,x_upper,x_fixed)
% SELECT_GRID_SAMPLES Select values for linearly spaced grid samples
%
%   X_GRID_SAMPLE = SELECT_GRID_SAMPLES(N_SAMPLES_PER_DIM,X_LOWER,X_UPPER) returns
%   a matrix X_GRID_SAMPLE where each row is a sample and the samples are
%   linearly spaced along all grid dimensions in the hyperrectangle defined by 
%   upper parameter values X_UPPER and lower parameter values X_LOWER, 
%   both row vectors.
%
%   X_GRID_SAMPLE = SELECT_GRID_SAMPLES(N_SAMPLES_PER_DIM,X_LOWER,X_UPPER,X_FIXED)
%   creates the same matrix of samples but only for variation along the
%   dimensions for which the parameter value in X_FIXED is NaN; for the 
%   other parameters, the values in the samples are fixed to the
%   corresponding value in X_FIXED.

% Number of parameters (i.e. dimensions in sample space)
n_params = size(x_lower,2);

if nargin < 4
    % By default, no parameters are fixed
    x_fixed = NaN(1,n_params);    
end

% True for dimensions for which parameter should have a fixed value
fix_dim = ~isnan(x_fixed);

% Indiced of fixed and variable dimensions
I_dim_fixed = find(fix_dim);
I_dim_var = find(~fix_dim);

% Number of fixed parameters
n_params_fixed = sum(fix_dim);

% Number of variable parameters
n_params_var = n_params - n_params_fixed;

% Spacing in linear grid for each dimension
dx = (x_upper - x_lower)./(n_samples_per_dim - 1);

% Set of parameter values in grid coordinates for each parameter
x_linspace = x_lower + (0:n_samples_per_dim-1)' * dx;

% Initialize samples  
combs = combinator(n_samples_per_dim,n_params_var,'p','r');
n_samples = size(combs,1);
x_grid_sample = zeros(n_samples,n_params);

% Set value of all parameters that are varied in the sample
for i_dim_var = 1:n_params_var
    i_dim = I_dim_var(i_dim_var);
    x_grid_sample(:,i_dim) = x_linspace(combs(:,i_dim_var),i_dim);
end

% Set value of all parameters that should be fixed
for i_dim_fixed = 1:n_params_fixed
    i_dim = I_dim_fixed(i_dim_fixed);
    x_grid_sample(:,i_dim) = x_fixed(i_dim);
end



