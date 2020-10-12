function [IC_alt_0,IC_alt_1,CB_0,CB_1,Rel_benefits] = prob_socioec_analysis(parameters_prob,parameters_fixed,socioeconomic_analysis_options)
% SOCIOECONOMIC_ANALYSIS Calculate socio-economic cost-benefit for HILP mitigation
%
% Carry out a socio-economic analysis of the contributions of interruption 
% costs due to HILP events and the costs of reducing the risk of HILP events.
%
% INPUTS:
%   parameters_prob:
%   Struct of parameters treated as having fixed values in this function (but that 
%   may have epistemic uncertainty)
%
%   parameters_fixed:
%   Struct of parameters treated as having fixed values in this function (but that 
%   may have epistemic uncertainty)
%  
%   socioeconomic_analysis_options:
%   Struct of options for the socio-economic analysis
%
% OUTPUTS:
%   IC_alt_0:
%   Column vector of samples of interruption costs for grid alternative 0
%
%   IC_alt_1:
%   Column vector of samples of interruption costs for grid alternative 1
%
%   CB_0:
%   Column vector of samples of cost minus benefits for grid alternative 0
%
%   CB_1:
%   Column vector of samples of cost minus benefits for grid alternative 1
%
%   Rel_benefits:
%   Vector of sampled relative cost-benefit of risk mitigating measures



%% Unpacking input parameters 

% Options
do_plot = socioeconomic_analysis_options.do_plot;
do_show_BC_and_not_CB = socioeconomic_analysis_options.do_show_BC_and_not_CB;
do_allow_max_1_event = socioeconomic_analysis_options.do_allow_max_1_event;
lambda_scaling = socioeconomic_analysis_options.lambda_scaling;
IC_scaling = socioeconomic_analysis_options.IC_scaling;
t_horizon = socioeconomic_analysis_options.t_horizon;
r_discount = socioeconomic_analysis_options.r_discount;
n_MC_sims = socioeconomic_analysis_options.n_MC_sims;
symmetric_uncertainty = socioeconomic_analysis_options.symmetric_uncertainty;

% Parameters treated as having fixed values in this function (but that may have epistemic uncertainty)
r_load_growth = parameters_fixed.r_load_growth;
r_incr_c = parameters_fixed.r_incr_c;
f_irred = parameters_fixed.f_irred;
l = parameters_fixed.l;

% Parameters treated probabilistically (with aleatory uncertainty)
r_input = parameters_prob.r_input;
c_input = parameters_prob.c_input;
P_input = parameters_prob.P_input;
C_0_input = parameters_prob.C_0_input;
C_1_input = parameters_prob.C_1_input;

%% Preprocessing

% Discount factor for interruption costs as a function of year
discount_factor = 1./( 1 + r_discount).^((1:t_horizon)'-1);

% Uncomment the following line to disable discounting
%discount_factor = ones(t_horizon,1);

% Factor accounting for expected increase in lost load and thus
% interruption costs due to load growth
load_growth_factor = (1 + r_load_growth).^(1:t_horizon)';

% Factor accounting for expected increase in valuation of lost load
c_growth_factor = (1 + r_incr_c).^(1:t_horizon)';

%% Sample grid costs (from triangular distribution)

if strcmp(C_0_input.model,'tri')
    C_0 = zeros(n_MC_sims,1);
    C_1 = zeros(n_MC_sims,1);
    for i_MC_sim = 1:n_MC_sims
        C_0(i_MC_sim) = rand_tri(C_0_input.lower,C_0_input.mode,C_0_input.upper);
        C_1(i_MC_sim) = rand_tri(C_1_input.lower,C_1_input.mode,C_1_input.upper);
    end
else
    C_0 = C_0_input.mode;
    C_1 = C_1_input.mode;
end

%% Sampling HILP events

% Logical arrays with value true if a HILP event happens for that year
% (row) in that MC iteration (column)
event_red_happens = false(n_MC_sims,t_horizon);
event_irred_happens = false(n_MC_sims,t_horizon);

% Sparse matrices with non-zero entries with value i_year for the elements
% (i_year, MC iteration) where a HILP event occurs 
year_event_red = sparse(t_horizon,n_MC_sims);
year_event_irred = sparse(t_horizon,n_MC_sims);
    
% Failure rate of HILP events (irreducible -- not eliminated by
% risk-mitigating measure)
lambda_irred = lambda_scaling * l * f_irred;

% Failure rate of HILP events (reducible -- is eliminated by
% risk-mitigating measure)
lambda_red = lambda_scaling * l * (1-f_irred);

% Annual probabilities of the occurrence of HILP events
prob = 1 - exp(- (lambda_red + lambda_irred) );

for i_sim = 1:n_MC_sims
   
    % Simulating HILP events
    for i_year = 1:t_horizon
        event_happens = (rand < prob);
        if event_happens
            event_is_irred = (rand < f_irred);
            if event_is_irred
                % The risk of this HILP event cannot be reduced
                % by the risk-reducing measure
                event_irred_happens(i_sim,i_year) = true;
                year_event_irred(i_year,i_sim) = i_year;
            else
                % The risk of this HILP event can be reduced by
                % the risk-reducing measure
                event_red_happens(i_sim,i_year) = true;
                year_event_red(i_year,i_sim) = i_year;
            end
            if do_allow_max_1_event
                % If assuming that a HILP event will trigger a risk-reducing
                % measure to avoid that it happens again
                break;
            end
        end        
    end
end

% Number of events during analysis period for each simulation
n_events_irred = sum(event_irred_happens,2);
n_events_red = sum(event_red_happens,2);

% Number of events in total
sum_events_irred = sum(n_events_irred);
sum_events_red = sum(n_events_red);

% Transform year of the event from sparse N_MC x t_horizon matrix to n_events_* x 1 array
year_event_irred = full(year_event_irred(event_irred_happens'));
year_event_red = full(year_event_red(event_red_happens'));

% Sample HILP events (irreducible)
IC_irred = IC_scaling * sample_interruption_cost([],P_input,r_input,c_input,sum_events_irred,symmetric_uncertainty,false);

% Sample HILP events (reducible)
IC_red = IC_scaling * sample_interruption_cost([],P_input,r_input,c_input,sum_events_red,symmetric_uncertainty,false);


%% Calculate and plot number of events in period

% Total number of HILP events in analysis period
n_events = n_events_irred + n_events_red;

uniqcount_events = uniqcount(n_events);

if do_plot
    fig = figure;
    axes1 = axes('Parent',fig);
    h_events = bar(uniqcount_events(:,1),uniqcount_events(:,2));
    set(h_events,'EdgeColor','k','FaceColor',[0.5 0.5 0.5])
    set(axes1,'YScale','log');
    
    ylim([0.1 n_MC_sims*4])
    xlabel('Number of HILP events during analysis period')
    ylabel('Number of Monte Carlo samples')
end

%% Calculating interruption costs

% Calculate interruption costs for each simulation / analysis horizon
IC_irred_per_sim = zeros(n_MC_sims,1);
IC_red_per_sim = zeros(n_MC_sims,1);

% HILP events (irreducible)
i_event_start = 1;
for i_sim = 1:n_MC_sims  
    i_event_end = i_event_start + n_events_irred(i_sim) - 1;
    I_events = i_event_start:i_event_end;

    % Interruption cost contributions for MC sim. accounting for load
    % growth and increased value of lost load
    IC_irred_in_sim = IC_irred(I_events) .* load_growth_factor(year_event_irred(I_events)) .* c_growth_factor(year_event_irred(I_events));
    
    % Interruption cost contributions for MC sim. accounting for discounting
    IC_irred_disc_in_sim = IC_irred_in_sim .* discount_factor(year_event_irred(I_events));    
    
    IC_irred_per_sim(i_sim) = sum(IC_irred_disc_in_sim);
    i_event_start = i_event_end + 1;      
end

% HILP events (reducible)
i_event_start = 1;
for i_sim = 1:n_MC_sims  
    i_event_end = i_event_start + n_events_red(i_sim) - 1;
    I_events = i_event_start:i_event_end;

    % Interruption cost contributions for MC sim. accounting for load
    % growth and increased value of lost load
    IC_red_in_sim = IC_red(I_events) .* load_growth_factor(year_event_red(I_events)) .* c_growth_factor(year_event_red(I_events));
    
    % Interruption cost contributions for MC sim. accounting for discounting
    IC_red_disc_in_sim = IC_red_in_sim .* discount_factor(year_event_red(I_events));    
    
    IC_red_per_sim(i_sim) = sum(IC_red_disc_in_sim);
    i_event_start = i_event_end + 1;      
end


%% Interruption costs for both grid development alternatives

IC_alt_0 = IC_irred_per_sim + IC_red_per_sim;
IC_alt_1 = IC_irred_per_sim;

%% Socio-economic cost-benefit for each alternative

if do_show_BC_and_not_CB
    % Impacts (benefits-costs) for zero alternative 
    CB_0 = -(C_0 + IC_irred_per_sim + IC_red_per_sim);
    
    % Impacts (benefits-costs) for risk-reducing measures 
    CB_1 = -(C_1 + IC_irred_per_sim);
else 
    % Socio-economic costs for zero alternative 
    CB_0 = C_0 + IC_irred_per_sim + IC_red_per_sim;
    
    % Socio-economic costs for risk-reducing measures 
    CB_1 = C_1 + IC_irred_per_sim;
end

%% Socio-economic impacts (benefits-costs) relative to zero alternative

if do_show_BC_and_not_CB
    Rel_benefits = CB_1 - CB_0;
else
    Rel_benefits = -(CB_1 - CB_0);
end
    
