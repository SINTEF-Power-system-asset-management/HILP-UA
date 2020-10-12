function [P_avg,t_start] = rand_P_time_series(P_time_series,r)  
%   RAND_P_TIME_SERIES Samples from load time series
%
%   P_AVG = RAND_P_TIME_SERIES(P_TIME_SERIES,R) returns the average load P_AVG
%   during an interruption event with duration R sampled from a load time 
%   series P_TIME_SERIES.
%
%   [P_AVG, T_START] = RAND_P_TIME_SERIES(P_TIME_SERIES,R) in addition
%   returns the time step T_START that the sample (interruption events)
%   starts.

    % Length of load time series (not doing any checks here)
    T = size(P_time_series,1);
    
    % Start time of interruption event
    t_start = ceil(rand*T);
    
    % End time of interruption event
    t_end = t_start + r - 1 ;
    
    % If end time is after end of time series, handle the time series as if
    % it was cyclic (NB: This assumes that the time series is an integer
    % number of years, but this is not checked now...)
    if t_end > T
        t_excess = t_end - T;
        I_t_interr = [t_start:T, 1:t_excess];
    else
        I_t_interr = t_start:t_end;
    end    
    
    % Calculate average load during interruption event
    P_time_series_interr = P_time_series(I_t_interr);
    P_avg = mean(P_time_series_interr);
end