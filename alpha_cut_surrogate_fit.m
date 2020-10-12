function [Z_fit_min,Z_fit_max] = alpha_cut_surrogate_fit(x_design,Z_samples,ub,lb,x_0,use_response_surface)
% ALPHA_CUT_SURROGATE_FIT Finds alpha-cut estimate based on surrogate model fit
%
%   [Z_FIT_MIN,Z_FIT_MAX] = ALPHA_CUT_SURROGATE_FIT(X_DESIGN,Z_SAMPLES,...
%       UB,LB,X_0,USE_RESPONSE_SURFACE)
%   returns estimate of the alpha cut for Z based on a response surface 
%   (surrogate) model fit to the function for Z using constrained
%   optimization within hyperrectangle. The response surface is a quadratic
%   multiple regression model that is defined by 15 coefficiencts and is
%   fitted based on n_params samples from the hyper-rectangle in the
%   epistemic uncertainty space. This space is defined by n_params_poss
%   parameters with epistemic uncertainty.
%
%
% INPUTS:
%
%   X_DESIGN: 
%   Design matrix for quadratic multiple regression (response surface
%   model); n_samples x 15 double array.
%
%   Z_SAMPLES:
%   Samples of Z from the hyperrectangle in the epistemic uncertainty space
%   that is used to fit the surrogate model; n_samples x 1 double array
%
%   UB:
%   Vertex defining the upper bounds of the hyper-rectangle; 
%   1 x 4 double array
%
%   LB:
%   Vertex defining the upper bounds of the hyper-rectangle; 
%   1 x 4 double array
%
%   X_0:
%   Initial point in the epistemic uncertainty space for the optimization 
%   solver; 1 x 4 double array.
%
%   USE_RESPONSE_SURFACE: 
%   Logical parameter that is false if response surface model is not to be 
%   used after all (optional; default: true)
%
%
% OUTPUTS:
%   Z_FIT_MIN: 
%   Lowest possible value of Z in the hyperrectangle
%
%   Z_FIT_MIN: 
%   Highest possible value of Z in the hyperrectangle 


if nargin < 6
    use_response_surface = true;
end

% Set up for optimization within hyperrectangle (fmincon)
A = [];
b = [];
Aeq = [];
beq = [];

if ~use_response_surface
    % If not using surrogate modelling (e.g. because there is not enough
    % samples to make a sensible fit)
    Z_fit_min = NaN;
    Z_fit_max = NaN;
else
    % Parameters of quadratic response surface model
    beta = (x_design' * x_design)\ (x_design' * Z_samples);
    
    % Defining functions to minimize to find maximum and minimum value of 
    % the response surface, respectively
    fun_max = @(x) -x2fx(x,'quadratic')*beta;
    fun_min = @(x) x2fx(x,'quadratic')*beta;
    
    % Find coordinates in epistemic uncertainty space for maximum and
    % minimum, respectively
    x_min = fmincon(fun_min,x_0,A,b,Aeq,beq,lb,ub);
    x_max = fmincon(fun_max,x_0,A,b,Aeq,beq,lb,ub);
    
    % Find maximum and minimum value of Z as fitted by the response surface
    % model
    Z_fit_min = x2fx(x_min,'quadratic')*beta;
    Z_fit_max = x2fx(x_max,'quadratic')*beta;
end