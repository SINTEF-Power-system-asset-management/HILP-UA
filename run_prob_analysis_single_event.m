% Runscript for probabilistic uncertainty analysis for the
% consequences of a single HILP event (given that it occurs)

%% Load default input parameters and settings

[parameters_prob, parameters_poss] = input_params_case_study;

c_input_default = parameters_prob.c_input;
l_input = parameters_poss.l_input;
P_input = parameters_prob.P_input;
r_input = parameters_prob.r_input;
symmetric_uncertainty = parameters_prob.symmetric_uncertainty;

%% Set parameters of the analysis 

% True if multiplying impact of HILP event in calculating interruption costs 
% (consequences) of HILP events (if that were the case, the analysis would 
% not be considering the consequences of a single HILP event)
do_consider_prob = false;

% Number of Monte Carlo samples
n_MC = 100000;

% Level for Value-at-Risk
alpha = 0.01;

% True if plotting Conditional Value-at-Risk (CVaR)
do_plot_CVaR = false;

% Label of x axis (consequence) in plots 
str_xlabel = 'Interruption costs due to a HILP event (NOK)';

%% Sample HILP events with time-dependent correlations (default settings)

c_input = c_input_default;

[IC_time_dep, ~, r_time_dep, P_time_dep, c_time_dep] = sample_interruption_cost(l_input,P_input,r_input,c_input,n_MC,symmetric_uncertainty,do_consider_prob);

% Store results as histogram for interruption costs
figure
h_time_dep = histogram(IC_time_dep,100);
format_histogram(h_time_dep)
xlabel(str_xlabel)
ylabel('Probability density function')

% Store marginal probability density function of specific interruption costs 
% to be able to sample from it without considering correlation with other
% input parameters
h_c_custom = histogram(c_time_dep,100,'Normalization','cdf');
xlabel('Specific interruption costs c (NOK/MWh)')
ylabel('Probability density function')

%% Sample HILP events without time-dependent correlations 

% Use sampled specific interruption costs as a custom probability distribution 
% instead of modelling time-dependence of specific interruption costs
% directly
c_input = c_input_default;
c_input.model = 'custom_distr';
c_input.custom_distr = h_c_custom;

[IC_no_corr,~, ~, P_no_corr, c_no_corr] = sample_interruption_cost(l_input,P_input,r_input,c_input,n_MC,symmetric_uncertainty,do_consider_prob);

% Store results as histogram for interruption costs
h_no_corr = histogram(IC_no_corr,100);
format_histogram(h_no_corr)
xlabel(str_xlabel)
ylabel('Probability density function')

%% Sample HILP events without uncertainty in specific interruption costs 

% Use "trivial" probability distribution for specific interruption cost 
% where the probability is 1 for getting the mean
c_input = c_input_default;
c_input.model = 'tri';
c_input.c_0 = mean(c_time_dep);
c_input.c_lower = c_input.c_0;
c_input.c_upper = c_input.c_0;

[IC_no_uncert,~, ~, P_no_uncert, c_no_uncert] = sample_interruption_cost(l_input,P_input,r_input,c_input,n_MC,symmetric_uncertainty,do_consider_prob);

% Store results as histogram for interruption costs
h_no_uncert = histogram(IC_no_uncert,100);
format_histogram(h_no_uncert)
xlabel(str_xlabel)
ylabel('Probability density function')

%% Plot PDF  with time-dependent correlations (default, benchmark)

plot_prob_distr(IC_time_dep,false,true,alpha,do_plot_CVaR)

%% Compare CDFs with and without time-dependent correlations

do_plot_log = false;
str_legend = {'With time-dependent correlations between $P$ and $c$','Without time-dependent correlations'};
compare_cdfs({h_time_dep,h_no_corr},alpha,do_plot_log,str_xlabel,str_legend,do_plot_CVaR)

%% Compare all three CDFs in a log-log plot

do_plot_log = true;
str_legend = {'With time-dependent correlations between $P$ and $c$','Without time-dependent correlations','Without uncertainty in c'};
compare_cdfs({h_time_dep,h_no_corr,h_no_uncert},alpha,do_plot_log,str_xlabel,str_legend,do_plot_CVaR)
