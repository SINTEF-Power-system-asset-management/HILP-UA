function socioeconomic_analysis_options = set_socioec_options()
% SET_SOCIOEC_OPTIONS Sets options to use for socio-economic analysis

% Analysis horizon (years)
t_horizon = 40;

% Discount rate
r_discount = 0.04;

% Number of Monte Carlo iterations of analysis (simulation) horizons
n_MC_sims = 100000;

% Scaling factor with which to multiply the failure rates of HILP events 
% -- to be used to check if approximate break-even analysis gives the
% expected results (i.e. break-even expected value) also when considering
% the distribution functions using the MC analysis below.
lambda_scaling = 1;

% Scaling factor for consequences (to be able to consider non-HILP event
% for comparison by setting it very low)
IC_scaling = 1;

% True if allowing at most 1 HILP event to occur (assuming risk-reducing 
% measures would have to be implemented if so happens, but assuming that the 
% cost of these measures is negligible
do_allow_max_1_event = true;

% True if plotting results
do_plot = false;

% True if defining socioeconomic impacts such that costs are shown as
% negative values, i.e. showing benefits-costs instead of costs-benefits
do_show_BC_and_not_CB = false;

% True if assuming symmetric probability distribution for uncertainties
% represented as triangular distributions
symmetric_uncertainty = false;

% Options for the probabilistic socio-economic analysis
socioeconomic_analysis_options.do_plot = do_plot;
socioeconomic_analysis_options.do_show_BC_and_not_CB = do_show_BC_and_not_CB;
socioeconomic_analysis_options.do_allow_max_1_event = do_allow_max_1_event;
socioeconomic_analysis_options.lambda_scaling = lambda_scaling;
socioeconomic_analysis_options.IC_scaling = IC_scaling;
socioeconomic_analysis_options.t_horizon = t_horizon;
socioeconomic_analysis_options.r_discount = r_discount;
socioeconomic_analysis_options.n_MC_sims = n_MC_sims;
socioeconomic_analysis_options.symmetric_uncertainty = symmetric_uncertainty;