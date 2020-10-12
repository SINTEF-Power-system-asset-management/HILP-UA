%% Create synthetic load time series 
% Create synthetic load time series to be able to test run the code and get
% sensible results without the real, measured and confidential load time
% series used in the case study. This code can only be run if the real load 
% time series is available and can be loaded to the workspace in P_input by 
% running input_params_case_study using use_real_load_data = true

% Real, measured load time series
P_real = P_input.P_time_series;

% Number of time steps (hours) of real load time series
n_time_steps = size(P_real,1);

% Timestamps of real load time series as a n_time_steps x 6 array on the
% MATLAB datevec format
datavec_P_time_series = P_input.datavec_P_time_series;

% Time vector (time resolution: 1 hour)
t = (1:n_time_steps)';

% Load variation parameters (MW) in simplistic load model
Delta_P_annual = 90;
Delta_P_diurnal = 52;
P_mean = 246.7;

% Time dependence parameters in simplistic load model
days_shift_annual = 31;
hours_shift_diurnal = 31;

% Extremely simplistic load model with purely sinusoidal annual variations
% and purely sinusoidal diurnal variations. This load variation reproduces
% the main characteristics (mean and standard deviation) of the real load 
% time series but does not fully reproduce the tails of the distribution
% function (i.e. the extreme variations)
P_synth = P_mean + cos((t-24*days_shift_annual)/8760*2*pi)*Delta_P_annual + ...
    cos((t-hours_shift_diurnal)/24*2*pi)*Delta_P_diurnal;

%% Write synthetic load time series to file
dlmwrite('case_data\P_time_series_synth.csv',[datavec_P_time_series, P_synth]);

%% Compare real and synthetic load time series

figure
hold on
plot(t,P_real,'k','LineWidth',0.5)
plot(t,P_synth,'r','LineWidth',0.5)
legend('Real, measured load time series','Synthetic load time series')
xlabel('Time (hour)')
ylabel('Load on transmission line (MW)')

mean_P_real = mean(P_real)
mean_P_synth = mean(P_synth)

std_P_real = std(P_real)
std_P_synth = std(P_synth)