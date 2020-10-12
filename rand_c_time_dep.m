function c = rand_c_time_dep(date_vec,r,c_ref_vec,f_c_ENS_per_month) 
%   RAND_C_TIME_DEP Samples specific interruption cost dependent on time
%
%   C = RAND_C_TIME_DEP(DATE_VEC,R,F_C_ENS_PER_MONTH) returns the 
%   specific interruption cost (NOK/MWh) for time step in a load time series 
%   DATE_VEC (on the format returned by the DATEVEC function), given the 
%   interruption duration R (in hours). C_REF_VEC is a column vector with 8760
%   elements where C_REF_VEC(R) gives the specific interruption cost in NOK/kWh 
%   for an interruption at the reference time with duration R hours. 
%   F_C_ENS_PER_MONTH is a column vector with 12 elements where each gives 
%   the average time-dependent correction factor for the specific interruption 
%   cost for the given month.

    if r > 8760
        error('Interruption durations > 8760 hours not supported')        
    end

    % Specific interruption cost at reference time for this duration
    c_ref = c_ref_vec( round(r) );
    
    % Find time-dependent correction factor for time (only information 
    % about the month in the current version)        
    m = date_vec(1,2);
    f_c_ENS = f_c_ENS_per_month(m);

    % Calculate specific interruption cost in NOK/MWh given the time and 
    % duration of the interruption event
    c = c_ref * f_c_ENS * 1000;
end