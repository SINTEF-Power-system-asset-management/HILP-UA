function x_lh_sample = select_lh_samples(n_samples,x_lower,x_upper)
% SELECT_LH_SAMPLES Select values for Latin Hypercube Samples
% 
%   X_LH_SAMPLE = SELECT_LH_SAMPLES(N_SAMPLES,X_LOWER,X_UPPER) returns
%   a matrix X_LH_SAMPLE with N_SAMPLES rows where each row contains
%   parameter values a Latin Hypercube Sample; the samples are uniformly 
%   distributed in the hyperrectangle defined by upper parameter values X_UPPER 
%   and lower parameter values X_LOWER, both row vectors.


% Number of parameters (i.e. dimensions in sample space)
n_params = size(x_lower,2);

% Uniformly distributed Latin Hypercube samples in a normalized coordinates
% (i.e. a hypercube with lower parameter values 0 and upper parameter
% values 1)
z = lhsdesign(n_samples,n_params);

% Rescale parameter values from normalized hypercube to the hypercube 
% defined in the inputs 
x_lh_sample = x_lower + z.*(x_upper - x_lower);
