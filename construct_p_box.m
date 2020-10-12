function cdf_values_p_box = construct_p_box(cdf_values,use_surrogate,x_poss_sample,I_sample_p_box)
% CONSTRUCT_P_BOX Constructs p-box from cumulative distribution functions
%   
%   CDF_VALUES_P_BOX = CONSTRUCT_P_BOX(CDF_VALUES,USE_SURROGATE,X_POSS_SAMPLE,I_SAMPLE_P_BOX)
%   takes a N_BINS x N_SAMPLES array of cumulative distribution functions
%   (CDFs) CDF_VALUES and estimates probability bounds (a "p-box") defined
%   by a lower CDF CDF_VALUES_P_BOX(:,1) and an upper CDF
%   CDF_VALUES_P_BOX(:,2), both column vectors with N_BINS elements. If
%   USE_SURROGATE_MODEL = true, a quadratic response surface (surrogate)
%   model is used in the alpha-cut method; if USE_SURROGATE_MODEL = false,
%   the lowest and highest values from the sampled CDFs are used to
%   estimate the alpha-cuts. X_POSS_SAMPLE is a N_SAMPLES x N_PARAM_POSS
%   array of the values of each of the N_PARAM_POSS parameters for each
%   sample. I_SAMPLE_P_BOX = [i_sample_lower,i_sample_higher,i_sample_0] are 
%   the indices of the samples from the epistemic uncertainty space 
%   corresponding to the lowest possible, highest possible and best-guess
%   value.

% Number of histogram bin edges for which the value of the CDFs are defined
n_bin_edges = size(cdf_values,1);

% Lower and upper bounds for the CDFs in the p-box
cdf_p_box_lower = zeros(n_bin_edges,1);
cdf_p_box_upper = zeros(n_bin_edges,1);

% Design matrix for quadratic response surface (quadratic regression) model 
x_design = x2fx(x_poss_sample,'quadratic');

% Indices of sample from epistemic uncertainty space corresponding to
% lowest possible, highest possible and best-guess value
i_sample_lower = I_sample_p_box(1);
i_sample_upper = I_sample_p_box(2);
i_sample_0 = I_sample_p_box(3);

% Lower possible, upper possible, and best-guess value of uncertain parameters
lb = x_poss_sample(i_sample_lower,:);
ub = x_poss_sample(i_sample_upper,:);
x_0 = x_poss_sample(i_sample_0,:);

% Estimate upper and lower possible value of Z for each value of Y
for i_bin_edge = 1:n_bin_edges
    % Samples of Z for this value of Y
    Z_samples = cdf_values(i_bin_edge,:)';
    if use_surrogate
        % Use surrogate model to estimate upper and lower possible value of Z
        [Z_fit_min,Z_fit_max] = alpha_cut_surrogate_fit(x_design,Z_samples,ub,lb,x_0,true);
        cdf_p_box_lower(i_bin_edge) = min([Z_samples;Z_fit_min]);
        cdf_p_box_upper(i_bin_edge) = max([Z_samples;Z_fit_max]);
    else
        % Select upper and lower value of Z from sampled values
        cdf_p_box_lower(i_bin_edge) = min(Z_samples);
        cdf_p_box_upper(i_bin_edge) = max(Z_samples);
    end
end

% CDF values have to be in the interval [0,1]
cdf_p_box_lower(cdf_p_box_lower < 0) = 0;
cdf_p_box_upper(cdf_p_box_upper > 1) = 1;

% Combine lower, upper and best-guess CDF to a n_bin_edges x 3 array
cdf_values_p_box = [cdf_p_box_lower, cdf_p_box_upper, cdf_values(:,i_sample_0)];

