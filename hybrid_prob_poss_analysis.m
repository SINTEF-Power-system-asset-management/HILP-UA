function [results_poss_distr, results_p_box] = hybrid_prob_poss_analysis(poss_settings,parameters_prob,parameters_poss,socioeconomic_analysis_options)
% HYBRID_PROB_POSS_ANALYSIS Performs hybrid probabilistic-possibilistic uncertainty analysis
%
%   [RESULTS_POSS_DISTR, RESULTS_P_BOX] = HYBRID_PROB_POSS_ANALYSIS(...
%       POSS_SETTINGS,PARAMETERS_PROB,PARAMETERS_POSS,SOCIOECONOMIC_ANALYSIS_OPTIONS)
%   performs a hybrid probabilistic-possibilistic uncertainty analysis
%   applied to a socio-economic cost-benefit analysis of power system
%   development alternatives and returns results on possibility
%   distributions and p-boxes (probability bounds).
%
% INPUTS:
%   POSS_SETTINGS:
%   Settings for the outer loop (possibilistic) loop of the hybrid 
%   probabilistic-possibilistic uncertainty analysis; struct with the
%   following fields:
%   .n_cuts - number of alpha cuts (including alpha = 0 and alpha = 1)
%   .n_lh_samples - number of latin hypercube samples
%   .use_surrogate_p_box - true if using surrogate modelling option for 
%       calculating p-boxes
%   .use_surrogate_poss_distr - true if using surrogate modelling option 
%       for calculating possibility distributions
%
%   PARAMETERS_PROB:
%   Input data for parameters with probabilistic uncertainty
%   representation
%
%   PARAMETERS_POSS:
%   Input data for parameters with possibilistic uncertainty
%   representation
%
%   SOCIOECONOMIC_ANALYSIS_OPTIONS:
%   Settings and parameters for the socio-economic cost-benefit analysis


% Unpack settings for the (outer loop of the) hybrid probabilistic-
% possibilistic uncertainty analysis
n_cuts = poss_settings.n_cuts;
n_lh_samples = poss_settings.n_lh_samples;
use_surrogate_poss_distr = poss_settings.use_surrogate_poss_distr;
use_surrogate_p_box = poss_settings.use_surrogate_p_box;

% Unpack input data for parameters with possibilistic uncertainty
% representation
l_input = parameters_poss.l_input;
f_irred_input = parameters_poss.f_irred_input;
r_incr_c_input = parameters_poss.r_incr_c_input;
r_load_growth_input = parameters_poss.r_load_growth_input;
n_params_poss_all = parameters_poss.n_params_poss_all;



%% Figure out how many parameters are treated possibilistically

% True if parameter is treated possibilistically;
is_param_incl = false(1,n_params_poss_all);

names_param_poss = cell(0);
n_params_poss = 0;
if strcmp(l_input.model,'tri')
    is_param_incl(1) = true;
    n_params_poss = n_params_poss + 1;
    names_param_poss = [names_param_poss, 'l_input'];
end
if strcmp(f_irred_input.model,'tri')
    is_param_incl(2) = true;
    n_params_poss = n_params_poss + 1;
    names_param_poss = [names_param_poss, 'f_irred_input'];
end    
if strcmp(r_load_growth_input.model,'tri')
    is_param_incl(3) = true;
    n_params_poss = n_params_poss + 1;
    names_param_poss = [names_param_poss, 'r_load_growth_input'];
end    
if strcmp(r_incr_c_input.model,'tri')
    is_param_incl(4) = true;
    n_params_poss = n_params_poss + 1;
    names_param_poss = [names_param_poss, 'r_incr_c_input'];
end    

%% Initialize possibilistic sampling

% Number of points for which the possibility function is to be evaluated
% (assuming triangular distributions)
n_points = n_cuts*2-1;

% Interval in \alpha between each cut
d_alpha = 1/(n_cuts-1);

% \alpha cuts, excluding \alpha = 0 and \alpha = 1
alpha_vec = 0:d_alpha:1;

% Vector of alpha values as used in representing the possibility function 
alpha_vec_poss = [alpha_vec alpha_vec(end-1:-1:1)];

% Initializing values from the possibility distributions for all uncertain 
% input parameters (columns) for set of points from the alpha cuts (rows)
x_poss_alpha = zeros(n_points,n_params_poss);

% Looping through all uncertain parameters and all \alpha cuts
for i_param_poss = 1:n_params_poss
    for i_cut = 1:n_cuts
        alpha = alpha_vec(i_cut);             
        x_i_core = parameters_poss.([names_param_poss{i_param_poss}]).core;
        x_i_lower = parameters_poss.([names_param_poss{i_param_poss}]).lower;
        x_i_upper = parameters_poss.([names_param_poss{i_param_poss}]).upper;
        
        % Set the lower and upper parameter values for the alpha cut
        x_poss_alpha(i_cut,i_param_poss) = x_i_lower + alpha*(x_i_core - x_i_lower);
        x_poss_alpha(n_points-i_cut+1,i_param_poss) = x_i_upper - alpha*(x_i_upper - x_i_core);
    end    
end

% The best-guess values for uncertain parameters
x_best_guess = x_poss_alpha(n_cuts,:);

%% Create samples

I_cut_for_sample = [];
x_poss_sample = zeros(0,n_params_poss);

% Looping over alpha cuts to add vertex samples
for i_cut = 1:n_cuts
   
    % Indices for upper and lower parameter values for alpha cut
    i_point_lower = i_cut;
    i_point_upper = n_points - i_cut + 1;
        
    % Select hyperrectangle vertex values for samples
    x_lower = x_poss_alpha(i_point_lower,:);
    x_upper = x_poss_alpha(i_point_upper,:);
    x_vertices = select_vertex_samples(n_params_poss,x_lower,x_upper);
    n_samples_cut = size(x_vertices,1);
    
    x_poss_sample = [x_poss_sample;x_vertices];
    I_cut_for_sample = [I_cut_for_sample; ones(n_samples_cut,1)*i_cut];
end

% The last sample in the set or vertex samples corresponds to the best-guess value
i_sample_0 = size(I_cut_for_sample,1);

% The index of the samples corresponding to the highest and lowest possible
% values of the uncertain parameters
i_sample_lower = find(ismember(x_poss_sample,min(x_poss_sample),'rows'));
i_sample_upper = find(ismember(x_poss_sample,max(x_poss_sample),'rows'));

% Add Latin Hypercube Samples
if n_lh_samples > 0
    x_lower = x_poss_alpha(1,:);
    x_upper = x_poss_alpha(end,:);
    x_lh_sample = select_lh_samples(n_lh_samples,x_lower,x_upper);    
    x_poss_sample = [x_poss_sample;x_lh_sample]; 
end

% Total number of possibilitic samples
n_samples = size(x_poss_sample,1);

% Figure out which samples should be considered when constructing
% possibility distribution for each alpha-cut (if so, idx_sample_in_cut = 1)
idx_sample_in_cut = false(n_samples,n_cuts);
for i_cut = 1:n_cuts
    x_lower = x_poss_alpha(i_cut,:);
    x_upper = x_poss_alpha(n_points - i_cut + 1,:);
    idx_sample_in_cut(:,i_cut) = all(x_poss_sample >= x_lower & x_poss_sample <= x_upper,2);
end

%% Do socioeconomic analysis

% Initialize sampling of possibilistic parameters
CB_0_samples_cell = cell(n_samples,1);
CB_1_samples_cell = cell(n_samples,1);
Mean_CB_0_samples = zeros(n_samples,1);
Mean_CB_1_samples = zeros(n_samples,1);
Prob_positive_net_benefit_samples = zeros(n_samples,1);
Mean_net_benefits_samples = zeros(n_samples,1);
Net_benefits_samples_cell =  cell(n_samples,1);

% Looping over vertices and sample results
for i_sample = 1:n_samples
    disp(['sample ' num2str(i_sample) ' of ' num2str(n_samples)])
    
    % Setting values of uncertain parameters that are treated as having
    % fixed values within probabilistic socio-economic analysis
    if strcmp(l_input.model,'tri')
        i_param = find(strcmpi(names_param_poss,{'l_input'}));
        parameters_fixed.l = x_poss_sample(i_sample,i_param);
    else
        parameters_fixed.l = l_input.core;
    end
    if strcmp(f_irred_input.model,'tri')
        i_param = find(strcmpi(names_param_poss,{'f_irred_input'}));
        parameters_fixed.f_irred = x_poss_sample(i_sample,i_param);
    else
        parameters_fixed.f_irred = f_irred_input.core;
    end
    if strcmp(r_incr_c_input.model,'tri')
        i_param = find(strcmpi(names_param_poss,{'r_incr_c_input'}));
        parameters_fixed.r_incr_c = x_poss_sample(i_sample,i_param);
    else
        parameters_fixed.r_incr_c = r_incr_c_input.core;
    end
    if strcmp(r_load_growth_input.model,'tri')
        i_param = find(strcmpi(names_param_poss,{'r_load_growth_input'}));
        parameters_fixed.r_load_growth = x_poss_sample(i_sample,i_param);
    else
        parameters_fixed.r_load_growth = r_load_growth_input.core;
    end
    
    % Run probabilistic socio-economic analysis
    [~,~,CB_0,CB_1,Rel_benefits] = prob_socioec_analysis(parameters_prob,parameters_fixed,socioeconomic_analysis_options);
    
    % Collect results
    Mean_net_benefits_samples(i_sample) = mean(Rel_benefits);
    Prob_positive_net_benefit_samples(i_sample) = mean(Rel_benefits > 0);
    Mean_CB_0_samples(i_sample) = mean(CB_0);
    Mean_CB_1_samples(i_sample) = mean(CB_1);
    CB_0_samples_cell{i_sample} = CB_0;
    CB_1_samples_cell{i_sample} = CB_1;
    Net_benefits_samples_cell{i_sample} = Rel_benefits;
end
    
%% Construct possibility distributions for output variables

% Initialize output variable arrays
Mean_rel_benefits_vec = zeros(n_points,1);
Prob_positive_net_benefit_vec = zeros(n_points,1);
Mean_CB_0_vec = zeros(n_points,1);
Mean_CB_1_vec = zeros(n_points,1);

Mean_CB_0_vec_vertex = zeros(n_points,1);
Mean_CB_1_vec_vertex = zeros(n_points,1);
Mean_rel_benefits_vec_vertex = zeros(n_points,1);
Mean_CB_0_vec_interior = zeros(n_points,1);
Mean_CB_1_vec_interior = zeros(n_points,1);
Mean_rel_benefits_vec_interior = zeros(n_points,1);

% Best-guess value of parameters with epistemic uncertainties is initial
% value if applying optimization in alpha-cut method
x_0 = x_best_guess;

% Looping over alpha cuts
for i_cut = 1:n_cuts
    disp(['alpha cut ' num2str(i_cut) ' of ' num2str(n_cuts) ])
   
    % Indices for upper and lower parameter values for alpha cut
    i_point_lower = i_cut;
    i_point_upper = n_points - i_cut + 1;
    
    % Set bounds for optimization within hyperrectangle (fmincon)
    lb = x_poss_alpha(i_point_lower,:);
    ub = x_poss_alpha(i_point_upper,:);
    
    % Indices of the samples to be considered for this alpha-cut
    I_samples_for_cut = find(idx_sample_in_cut(:,i_cut));
    
    % Design matrix for quadratic multiple regression (response surface model)
    % (NB: needs statistics toolbox)
    x_design = x2fx(x_poss_sample(I_samples_for_cut,:),'quadratic');
    
    % True if using response surface fitting results (only if more samples
    % than the number of parameters for the response surface model)
    use_response_surface = ( numel(I_samples_for_cut) > size(x_design,2) ) && use_surrogate_poss_distr;
    
    % Results for expected net benefits of risk-reducing measures
    %  using response surface model
    Z_samples = Mean_net_benefits_samples(I_samples_for_cut);
    [Y_fit_min,Y_fit_max] = alpha_cut_surrogate_fit(x_design,Z_samples,ub,lb,x_0,use_response_surface);
    
    %  using sampled points in epistemic uncertainty space   
    [~,i_sample_for_cut_lower] = min(Mean_net_benefits_samples(I_samples_for_cut));
    [~,i_sample_for_cut_upper] = max(Mean_net_benefits_samples(I_samples_for_cut));
    i_sample_lower = I_samples_for_cut(i_sample_for_cut_lower);
    i_sample_upper = I_samples_for_cut(i_sample_for_cut_upper);
    
    %  choosing min/max from response surface model and vertices and other samples)
    Mean_rel_benefits_vec(i_point_lower) = min( Mean_net_benefits_samples(i_sample_lower), Y_fit_min);  
    Mean_rel_benefits_vec(i_point_upper) = max( Mean_net_benefits_samples(i_sample_upper), Y_fit_max);

    Mean_rel_benefits_vec_vertex(i_point_lower) = Mean_net_benefits_samples(i_sample_lower);
    Mean_rel_benefits_vec_vertex(i_point_upper) = Mean_net_benefits_samples(i_sample_upper); 
    Mean_rel_benefits_vec_interior(i_point_lower) = Y_fit_min;
    Mean_rel_benefits_vec_interior(i_point_upper) = Y_fit_max;    
    
    
    % Results for probability of positive net benefits of risk-reducing measures
    %  using response surface model
    Z_samples = Prob_positive_net_benefit_samples(I_samples_for_cut);
    [Y_fit_min,Y_fit_max] = alpha_cut_surrogate_fit(x_design,Z_samples,ub,lb,x_0,use_response_surface);

    %  using sampled points in epistemic uncertainty space  
    [~,i_sample_for_cut_lower] = min(Prob_positive_net_benefit_samples(I_samples_for_cut));
    [~,i_sample_for_cut_upper] = max(Prob_positive_net_benefit_samples(I_samples_for_cut));
    i_sample_lower = I_samples_for_cut(i_sample_for_cut_lower);
    i_sample_upper = I_samples_for_cut(i_sample_for_cut_upper);
    %  choosing min/max from response surface model and vertices and other samples)    
    Prob_positive_net_benefit_vec(i_point_lower) = min( Prob_positive_net_benefit_samples( i_sample_lower), Y_fit_min);
    Prob_positive_net_benefit_vec(i_point_upper) = max( Prob_positive_net_benefit_samples( i_sample_upper), Y_fit_max); 
    
    
    % Results for expected cost-benefits of alternative 0
    %  using response surface model
    Z_samples = Mean_CB_0_samples(I_samples_for_cut);
    [Y_fit_min,Y_fit_max] = alpha_cut_surrogate_fit(x_design,Z_samples,ub,lb,x_0,use_response_surface);

    %  using sampled points in epistemic uncertainty space      
    [~,i_sample_for_cut_lower] = min(Mean_CB_0_samples(I_samples_for_cut));
    [~,i_sample_for_cut_upper] = max(Mean_CB_0_samples(I_samples_for_cut));
    i_sample_lower = I_samples_for_cut(i_sample_for_cut_lower);
    i_sample_upper = I_samples_for_cut(i_sample_for_cut_upper); 
    %  choosing min/max from response surface model and vertices and other samples)
    Mean_CB_0_vec(i_point_lower) = min( Mean_CB_0_samples(i_sample_lower), Y_fit_min);
    Mean_CB_0_vec(i_point_upper) = max( Mean_CB_0_samples(i_sample_upper), Y_fit_max);
    
    Mean_CB_0_vec_vertex(i_point_lower) = Mean_CB_0_samples(i_sample_lower);
    Mean_CB_0_vec_vertex(i_point_upper) = Mean_CB_0_samples(i_sample_upper);    
    Mean_CB_0_vec_interior(i_point_lower) = Y_fit_min;
    Mean_CB_0_vec_interior(i_point_upper) = Y_fit_max;
  
    
    % Results for expected cost-benefits of alternative 1
    % using response surface model
    Z_samples = Mean_CB_1_samples(I_samples_for_cut);
    [Y_fit_min,Y_fit_max] = alpha_cut_surrogate_fit(x_design,Z_samples,ub,lb,x_0,use_response_surface);        

    %  using sampled points in epistemic uncertainty space  
    [~,i_sample_for_cut_lower] = min(Mean_CB_1_samples(I_samples_for_cut));
    [~,i_sample_for_cut_upper] = max(Mean_CB_1_samples(I_samples_for_cut));
    i_sample_lower = I_samples_for_cut(i_sample_for_cut_lower);
    i_sample_upper = I_samples_for_cut(i_sample_for_cut_upper);
    %  choosing min/max from response surface model and vertices and other samples)
    Mean_CB_1_vec(i_point_lower) = min( Mean_CB_1_samples(i_sample_lower), Y_fit_min);
    Mean_CB_1_vec(i_point_upper) = max( Mean_CB_1_samples(i_sample_upper), Y_fit_max);
    
    Mean_CB_1_vec_vertex(i_point_lower) = Mean_CB_1_samples(i_sample_lower);
    Mean_CB_1_vec_vertex(i_point_upper) = Mean_CB_1_samples(i_sample_upper); 
    Mean_CB_1_vec_interior(i_point_lower) = Y_fit_min;
    Mean_CB_1_vec_interior(i_point_upper) = Y_fit_max;
end

%% Constructing p-boxes

% Samples of epistemic uncertainty space used to set the lower and upper
% boundaries as well as the initial value of the bounded optimization in
% the alpha-cut method
I_samples_p_box = [i_sample_lower, i_sample_upper, i_sample_0];

% Net socio-economic benefits of risk-reducing measures
x_lims = [-2E9 16E9];
[cdf_values,bin_edges_Net_benefits] = interp_cdf(Net_benefits_samples_cell,x_lims);
cdf_values_Net_benefits = construct_p_box(cdf_values,use_surrogate_p_box,x_poss_sample,I_samples_p_box);

% Total socio-economic costs for zero alternative
[cdf_values,bin_edges_CB] = interp_cdf(CB_0_samples_cell,x_lims);
cdf_values_CB_0 = construct_p_box(cdf_values,use_surrogate_p_box,x_poss_sample,I_samples_p_box);

% Total socio-economic costs for the risk-reducing measures
[cdf_values,~] = interp_cdf(CB_1_samples_cell,x_lims);
cdf_values_CB_1 = construct_p_box(cdf_values,use_surrogate_p_box,x_poss_sample,I_samples_p_box);

%% Packaging results

% Unpacking results for possibility distributions
results_poss_distr.alpha_vec_poss = alpha_vec_poss;
results_poss_distr.Mean_rel_benefits_vec = Mean_rel_benefits_vec;
results_poss_distr.Mean_rel_benefits_vec_vertex = Mean_rel_benefits_vec_vertex;
results_poss_distr.Mean_CB_0_vec = Mean_CB_0_vec;
results_poss_distr.Mean_CB_0_vec_vertex = Mean_CB_0_vec_vertex;
results_poss_distr.Mean_CB_1_vec = Mean_CB_1_vec;
results_poss_distr.Mean_CB_1_vec_vertex = Mean_CB_1_vec_vertex;
results_poss_distr.Prob_positive_net_benefit_vec = Prob_positive_net_benefit_vec;

% Unpacking results for p-boxes
results_p_box.bin_edges_Mean_rel_benefits = bin_edges_Net_benefits;
results_p_box.bin_edges_CB = bin_edges_CB;
results_p_box.cdf_values_CB_0 = cdf_values_CB_0;
results_p_box.cdf_values_CB_1 = cdf_values_CB_1;
results_p_box.cdf_values_Net_benefits = cdf_values_Net_benefits;
