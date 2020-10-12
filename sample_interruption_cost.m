function [IC, l, r, P, c] = sample_interruption_cost(l_input,P_input,r_input,c_input,n_samples,symmetric_uncertainty,do_consider_prob)
% SAMPLE_INTERRUPTION_COST Make a number of samples of interruption costs
%
% INPUTS:
%   l_input:
%   Input parameters parameterizing distribution function for failure rate
%   for this type of HILP events. (Not used if do_consider_prob = false,
%   and may in that case e.g. be set to be empty.)
%
%   P_input:
%   Input parameters parameterizing distribution function for power lost
%   for this type of HILP events
%
%   r_input:
%   Input parameters parameterizing distribution function for interruption
%   duration for this type of HILP events
%
%   c_input:
%   Input parameters parameterizing distribution function for specific
%   interruption costs for this type of HILP events
%
%   n_samples:
%   Number of Monte Carlo samples
%
%   symmetric_uncertainty:
%   True if assuming symmetric probability distribution for uncertainties
%   represented as triangular distributions
%
%   do_consider_prob
%   True if considering the probability in the calculation, for estimating
%   uncertainty in and calculating the expected interruption cost; false 
%   if calculating probability distribution of consequence (in which case 
%   uncertainties are aleatory)
%
% OUTPUTS:
%   IC:
%   n_samples x 1 array of sampled interruption cost values (NOK)
%
%   l: 
%   n_samples x 1 array of sampled \lambda values (1/year)
%
%   r
%   n_samples x 1 array of sampled interruption duration values (h)
%
%   P
%   n_samples x 1 array of sampled average interrupted power values (MW)
%
%   c:
%   n_samples x 1 array of sampled specific interruption cost values (NOK/MWh)


% If time series input for load, make sure one does not set up and use the
% distribution function parameters
if strcmp(P_input.model,'time_series')
    set_up_P_distr = false;
    P_time_series = P_input.P_time_series;
else
    set_up_P_distr = true;
end

% If using time-dependent specific interruption costs, make sure one does 
% not set up and use the distribution function parameters.
if strcmp(c_input.model,'time_dep')
    if ~strcmp(P_input.model,'time_series')
        % (If also using time series input for load, this captures correlations; 
        % if not, time-dependent interruption costs make little sense...) 
        error('Time-dependent specific interruption costs only supported for time-series input for the load')
    end
    set_up_c_distr = false;
    c_ref_vec = c_input.c_ref_vec; 
    f_c_ENS_per_month = c_input.f_c_ENS_per_month;
elseif strcmp(c_input.model,'custom_distr')    
    custom_distr = c_input.custom_distr;
    set_up_c_distr = false;
else
    set_up_c_distr = true;
end

if do_consider_prob
    l_0 = l_input.core;
end
if set_up_P_distr
    P_0 = P_input.P_0;
end
r_0 = r_input.r_0;
if set_up_c_distr
    c_0 = c_input.c_0;
end

if symmetric_uncertainty
    if do_consider_prob
        l_upper = l_0 + l_input.d_l;
        l_lower = l_0 - l_input.d_l;
    end

    if set_up_P_distr
        P_upper = P_0 + P_input.d_P;
        P_lower = P_0 - P_input.d_P;
    end
    
    r_upper = r_0 + r_input.d_r;
    r_lower = r_0 - r_input.d_r;    
    
        
    if set_up_c_distr
        c_upper = c_0 + c_input.d_c;
        c_lower = c_0 - c_input.d_c;
    end
else
    if do_consider_prob
        l_upper = l_input.l_upper;
        l_lower = l_input.l_lower;
    end

    if set_up_P_distr
        P_upper = P_input.P_upper;
        P_lower = P_input.P_lower;
    end
    
    r_upper = r_input.r_upper;
    r_lower = r_input.r_lower;
    
    if set_up_c_distr
        c_upper = c_input.c_upper;
        c_lower = c_input.c_lower;
    end
end

if isfield(c_input,'c_r')
    c_r = c_input.c_r;
else
    c_r = 0;
end


l = zeros(n_samples,1);
r = zeros(n_samples,1);
P = zeros(n_samples,1);
c = zeros(n_samples,1);
IC = zeros(n_samples,1);

% Sample from selected probability distributions
for i = 1:n_samples
    if do_consider_prob
        switch l_input.model
            case 'normal'
                l(i) = max(0, l_0 + d_l * randn);
            case 'tri'
                l(i) = rand_tri(l_lower,l_0,l_upper);
            otherwise
                error(['Unsupported distribution: ' model_l])
        end
    else
        l(i) = 1;
    end
     switch r_input.model
        case 'normal'
            r(i) = max(0, r_0 + d_r * randn);
        case 'tri'
            r(i) = rand_tri(r_lower,r_0,r_upper);
        case 'tri_corr'
            r(i) = rand_tri_corr(r_0,d_r,P(i),P_0,d_P,rho_r_P);
        otherwise
            error(['Unsupported distribution: ' model_r])
    end       
    switch P_input.model
        case 'normal'
            P(i) = max(0, P_0 + d_P * randn);
        case 'tri'
            P(i) = rand_tri(P_lower,P_0,P_upper);
        case 'time_series'
            [P(i), t] = rand_P_time_series(P_time_series,r(i));
        otherwise
            error(['Unsupported distribution: ' model_P])
    end
    switch c_input.model
        case 'normal'
            c(i) = max(0, c_0 + d_c * randn);
        case 'tri'
            c(i) = rand_tri(c_lower,c_0,c_upper);
        case 'time_dep'
            date_vec_t = P_input.datavec_P_time_series(t,:);
            c(i) = rand_c_time_dep(date_vec_t,r(i),c_ref_vec,f_c_ENS_per_month);       
        case 'custom_distr'
            c(i) = rand_custom_distr(custom_distr);
        otherwise
            error(['Unsupported distribution: ' c_input.model])            
    end
    IC(i) = interr_cost_model(l(i),P(i),r(i),c(i),c_r);
end
disp('Done sampling interruption cost')