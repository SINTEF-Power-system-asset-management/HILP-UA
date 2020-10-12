function [parameters_prob, parameters_poss] = input_params_case_study
% INPUT_PARAMS_CASE_STUDY Returns input parameter values for case study
%
%   [PARAMETERS_PROB, PARAMETERS_POSS] = INPUT_PARAMS_CASE_STUDY returns
%   input data for parameters with probabilistic (PARAMETERS_PROB) and 
%   possibilistic (PARAMETERS_POSS) uncertainty representation for case
%   study on accounting for uncertainties due to HILP events in power
%   system development.
%
% See also: SET_SOCIOEC_OPTIONS.M

% Probability model for parameters with probabilistic representation; the
% default for the case study is sampling from a time series for the
% interrupted power P, sampling from a triangular distribution for the
% interruption duration r and using time-dependent factors to model the
% samples of specific interruption costs c (which will then correlate with P).
P_input.model = 'time_series';
r_input.model = 'tri';
c_input.model = 'time_dep';

% Best-guess (or most likely) values when using triangular probability
% distributions
P_input.P_0 = 250;           % Interrupted power (MW)
r_input.r_0 = 7*24;          % Interruption duration (hours)
c_input.c_0 = 55*1E3;        % Specific interruption costs (NOK/MWh)

% Upper and lower limits of triangle when using triangular probability
% distributions
P_input.P_upper = 450;
P_input.P_lower = 100;
r_input.r_upper = 21*24;
r_input.r_lower = 3.5*24;
c_input.c_upper = 60E3;
c_input.c_lower = 40E3;


% If using time series input data for load (reset other input parameters)
if strcmp(P_input.model,'time_series')

    use_real_load_data = false;
    % True if reading real, measured load data (confidential and not made
    % publicly available), false if reading synthetic load time series 
    % simplisticly generated to resemble the main characteristics of the
    % measured data
    
    if use_real_load_data
        % File with time series load data, with range of Excel cells and column
        % with the load data to use
        filename_load_data = 'data.xlsx';
        file_range = 'C6:F112328';
        col_load = 4;
        
        % Time stamp of first and last entry in load data set
        % (note that last hour value is 16, but datevec uses the time this hour
        % starts, i.e. 15 instead of 16)
        datevec_first = [2007 1 1 0 0 0];
        datevec_last = [2019 10 25 15 0 0];
        
        % First and last (complete) year to extract from data set
        year_first = 2013;
        year_last = 2018;
        
        % Read load time series and corresponding time stamps from file
        [P_input.P_time_series, P_input.datavec_P_time_series] = read_load_data(filename_load_data,file_range,col_load,datevec_first,datevec_last,year_first,year_last);
    else
        filename_load_data = 'case_data\P_time_series_synth.csv';
        load_data = dlmread(filename_load_data);
        P_input.datavec_P_time_series = load_data(:,1:6);
        P_input.P_time_series = load_data(:,7);        
    end
    
    % Remove fields in uncertainty description that are not used for this
    % probability model
    P_input = rmfield(P_input,{'P_upper','P_lower','P_0'});
end

% If using time-dependent specific interruption costs (reset other input parameters)
if strcmp(c_input.model,'time_dep')
    % Configuration file for OPAL case set up for specifying interruption cost function
    inputfile_OPAL = fullfile(pwd,'case_data\hilp_uncertainty_case.toml');
    inputfile_c_ref_vec = fullfile(pwd,'case_data\c_ref_vec.csv');
    inputfile_f_c_ENS_per_month = fullfile(pwd,'case_data\f_c_ENS_per_month.csv');
    
    % True if using functionality in OPAL prototype implementation 
    % (MATLAB package developed by SINTEF Energy Research, not publicly 
    % available) to calculate interruption cost data according to the
    % Norwegian Cost of Energy Not Supplied (KILE) scheme
    use_OPAL = false;
    
    % Hard coding data for composition of customer types (at the one bus considered)
    cust_type_comp = [0.28 0.11 0.41 0.04 0.14 0.02];
    
    if use_OPAL
        % Set up interruption cost data etc. using functionality in OPAL
        % prototype implementation
        [c_input.c_ref_vec, c_input.f_c_ENS_per_month] = set_up_intcost_data(inputfile_OPAL,P_input,cust_type_comp);
    else
        % Load pre-processed case input files for interruption cost data
        c_input.c_ref_vec = dlmread(inputfile_c_ref_vec);
        c_input.f_c_ENS_per_month = dlmread(inputfile_f_c_ENS_per_month);
    end
    
    % Remove fields in uncertainty description that are not used for this
    % probability model
    c_input = rmfield(c_input,{'c_upper','c_lower','c_0'});
end

%% Grid development costs

% Grid investment cost data 
C_0 = 950 * 1E6;   % Grid cost of zero alternative
C_1 = 1150 * 1E6;   % Grid-cost of risk-reducing measures

% Best-guess (or most likely) values (i.e. the mode of the triangular distribution)
C_0_input.mode = C_0;
C_1_input.mode = C_1;

% Uncertainty interval
relative_uncertainty_C = 0.1; 
C_0_input.lower = C_0 * (1 - relative_uncertainty_C);
C_0_input.upper = C_0 * (1 + relative_uncertainty_C);
C_1_input.lower = C_1 * (1 - relative_uncertainty_C);
C_1_input.upper = C_1 * (1 + relative_uncertainty_C);

% Distribution function model ('tri' for triangular and 'none' for no
% probability distribution, i.e deterministic value)
C_0_input.model = 'tri';
C_1_input.model = 'tri';

%% Parameters with uncertainty that are represented possibilistically
% (Parameters associated with epistemic uncertainty; for outer loop of
% hybrid probabilistic-possibilistic uncertainty analysis)

% Frequency of occurrence of HILP events (rate per year)
l_input.core = 0.001;         
l_input.upper = 0.01;
l_input.lower = l_input.core / 10;
l_input.model = 'tri';

% Fraction of irreducible HILP events
f_irred_input.core = 0.1;
f_irred_input.lower = 0;
f_irred_input.upper = 0.5;
f_irred_input.model = 'tri';

% Percentage increase in specific interruption costs (valuation of security of supply) per year
r_incr_c_input.core = 0.013;
r_incr_c_input.lower = 0;
r_incr_c_input.upper = 0.02;
r_incr_c_input.model = 'tri';

% Expected percentage-wise load growth per year
r_load_growth_input.core = (1.22)^(1/23)-1;
r_load_growth_input.lower = (1.15)^(1/23)-1;
r_load_growth_input.upper = (1.3)^(1/23)-1;
r_load_growth_input.model = 'tri';

% Number of parameters that can be treated possibilistically
n_params_poss_all = 4;

%% Other parameters

% True if assuming symmetric probability distribution for uncertainties
% represented as triangular distributions
symmetric_uncertainty = false;

%% Package input parameters for probabilistic uncertainty analysis

parameters_prob = struct();
parameters_prob.r_input = r_input;
parameters_prob.P_input = P_input;
parameters_prob.c_input = c_input;
parameters_prob.C_0_input = C_0_input;
parameters_prob.C_1_input = C_1_input;
parameters_prob.symmetric_uncertainty = symmetric_uncertainty;

%% Package input parameters for possibilistic uncertainty analysis

parameters_poss = struct();
parameters_poss.l_input = l_input;
parameters_poss.f_irred_input = f_irred_input;
parameters_poss.r_incr_c_input = r_incr_c_input;
parameters_poss.r_load_growth_input = r_load_growth_input;
parameters_poss.n_params_poss_all = n_params_poss_all;

