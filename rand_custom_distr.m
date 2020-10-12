function y = rand_custom_distr(custom_distr)
% TRI_CORR Sample from triangular distribution
%
%   Y = RAND_CUSTOM_DISTR(CUSTOM_DISTR) returns a variable Y drawn 
%   from a cumulative distribution function (CDF) specified by the Values 
%   and Bin_Edges fields of the CUSTOM_DISTR struct (cf. also the HISTOGRAM
%   function).

% Extract CDF data; the CDF is defined as a function F(y), where the bin
% edges are along the y axis and the CDF values along the F axis
Values = custom_distr.Values;
BinEdges = custom_distr.BinEdges;

% Remove CDF values that are identical (i.e. combine bins)
[F_edges, idx_unique] = unique([0, Values],'last');

% Finding values for F and y used to find the interpolation function F(y) 
% and in turn the inverse y(F)
F_mid = ( F_edges(1:end-1) + F_edges(2:end) ) / 2;
y_mid = BinEdges(idx_unique(1:end-1));

% Draw random value of F
r = rand;
F_interp = r;

% Find corresponding value of y(F) from interpolation
y_interp = interp1(F_mid,y_mid,F_interp,'linear','extrap');
y = y_interp;
