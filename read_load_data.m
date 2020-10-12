function [P_time_series,datevec_P_time_series] = read_load_data(filename,range,col_load,datevec_first,datevec_last,year_first,year_last)
% READ_LOAD_DATA Reads load data from file
%
% INPUTS:
%   filename: 
%   Excel file with load time series
%
%   range:
%   Cell range in Excel with data; 1st column is hour of day (1-24)
%
%   col_load:
%   Column of the array extracted from range that contains the load time
%   series to use.
%
%   datevec_first:
%   Date vector (see function datevec) for first data point in range
%
%   datevec_first:
%   Date vector (see function datevec) for last data point in range
%
%   year_first:
%   First (full) year to extract from time series
%
%   year_last:
%   Last (full) year to extract from time series
%
% OUTPUTS:
%   P_time_series:
%   Column vector with hourly load time series in MW
%
%   datevec_P_time_series:
%   n_data x 6 array of date (row) vectors (see function datevec)
%   specifying time stamps for load time series entries


data = xlsread(filename,1,range);
n_cols = size(data,2);
n_data = size(data,1);
hour = data(:,1);

% Find gaps where an hour is missing (only gaps of a single hour was found
% for this time series)
idx_gap = [(diff(hour) == 2); false];
I_gap_orig = find(idx_gap);
n_gaps = numel(I_gap_orig);

% Create new data matrix with gaps with value zero
n_data_fill = n_data + n_gaps;
I_gap_fill = I_gap_orig + (1:n_gaps)';
[~,idx_gap_fill] = ismember((1:n_data_fill)',I_gap_fill);
data_fill = zeros(n_data_fill,n_cols);
data_fill(~idx_gap_fill,:) = data(:,:);

% Fill gaps by linear interpolation
x = find(~idx_gap_fill);
v = data_fill(~idx_gap_fill,:);
xq = (1:n_data_fill)';
data_fill = interp1(x,v,xq,'linear');

%% Parse .xlsx data and form matrices for the timestamp and for the load

% Corresponding time stamps 
date_vec = datevec(datenum(datevec_first):1/24:datenum(datevec_last));

% Extract data for selected years
t_start = find( date_vec(:,1) == year_first, 1, 'first');
t_end = find( date_vec(:,1) == year_last, 1, 'last');

% Define load on double circuit (for Excel data, gaps are already filled
% and we assume that peaks are already filtered out)
P_time_series = data_fill(t_start:t_end,col_load);

% Corresponding time stamps
datevec_P_time_series = datevec(datenum(date_vec(t_start,:)):1/24:datenum(date_vec(t_end,:)));
