function x_vertices = select_vertex_samples(n_params,x_lower,x_upper)
% SELECT_VERTEX_SAMPLES Select values for hyperrectangle vertex sample
%
%   X_VERTICES = SELECT_VERTEX_SAMPLES(N_PARAMS,X_LOWER,X_UPPER) returns
%   a matrix X_VERTICES with one column for each of the N_PARAMS parameters
%   one is sampling from where each row corresponds to a vertex of the
%   hyperrectangle defined by upper parameter values X_UPPER and lower
%   parameter value X_LOWER, each a N_PARAMS-element row vector.

% Initialize arrays and cell arrays for storing results for samples of
% epistemic uncertainty realizations
n_vertices = 2^n_params;

% Vertex combination matrix with samples along rows and parameters
% along columns; value 1 meaning lower-value vertex and 2 meaning
% upper-value vertex for that parameter and sample
combs = combinator(2,n_params,'p','r');

% Initialize samples of uncertain parameters from possibility distributions
x_vertices = zeros(n_vertices,n_params);
    
x_values = [x_lower; x_upper];

% Loop over uncertain parameters and over vertices
for i_param_poss = 1:n_params
    for i_poss_sample = 1:n_vertices
        % Which point on the possibility distribution does this vertex
        % corresponds to (either upper or lower value of alpha cut)
        i_point = combs(i_poss_sample,i_param_poss);
        
        % The upper/lower value of alpha cut for this parameter and sample
        x_vertices(i_poss_sample,i_param_poss) = x_values(i_point,i_param_poss);
    end
end

% Make sure only unique vertices are returned if the hyperrectangle is
% "degenerate" along one of its directions, e.g. if a cube turns out to be 
% only a square (2D) or a single point (1D)
x_vertices = unique(x_vertices,'rows');