function [cdf_values_interp,bin_edges_interp] = interp_cdf(results_cell,x_lims,n_bins)
% CDF_VALUES_INTERP Constructs a set of CDFs F(x) interpolated to the same x's
% 
%   [CDF_VALUES_INTERP,BIN_EDGES_INTERP] = INTERP_CDF(RESULTS_CELL,X_LIMS)
%   takes a set of sampled values in RESULTS_CELL, which is a N_CDFS x 1 
%   cell array of sample vectors, and constructs Cumulative Distributions 
%   Functions CDF_VALUES_INTERP (N_BINS-1 x N_CDFS array) interpolated so
%   that all CDFs are defined for the same histogram bin edges BIN_EDGES_INTERP
%   (N_BIN-1 x 1 array). The range of values for which the p-box will be 
%   calculated is given by the 1 x 2 array x_lims (optional). N_BINS is the
%   number of bins of the histograms used to construct empirical CDFs
%   (optional; default: 360).

% If argument x_lims is not provided, we will try find it from the
% sampled values in results_cell
if nargin < 2
    do_find_x_lims = true;
    x_lims = [Inf -Inf];
else
    do_find_x_lims = false;
end
    
% Number of bins in histograms used to form CDFs
if nargin < 3
    n_bins = 360;
end

n_cdfs = numel(results_cell);
bin_edges = cell(n_cdfs,1);
cdfs = cell(n_cdfs,1);

% Find highest and lowest values in results_cell
if do_find_x_lims
    for i_hist = 1:n_cdfs
        x_max_i = max(results_cell{i_hist});
        if x_max_i > x_lims(2)
            x_lims(2) = x_max_i;
        end        
        
        x_min_i = min(results_cell{i_hist});
        if x_min_i < x_lims(1)
            x_lims(1) = x_min_i;
        end
    end
end

% Find bin width for histograms used to estimate CDFs
bin_width = ( x_lims(2) - x_lims(1) )/(n_bins-1);
 
% Construct CDFs from histograms
for i_hist = 1:n_cdfs
    y = results_cell{i_hist};
    y(y > x_lims(2)) = x_lims(2);
    [cdfs{i_hist}, bin_edges{i_hist}] = histcounts(y,'Normalization','cdf','BinWidth',bin_width,'BinLimits',x_lims); 
end

% Number of CDFs to interpolate
n_cdfs = numel(cdfs);

% X values for interpolation of CDFs
bin_edges_interp = unique([bin_edges{:}]);
n_bins = numel(bin_edges_interp)-1;
bin_edges_interp = bin_edges_interp(1:n_bins)';
cdf_values_interp = zeros(n_bins,n_cdfs);

% Interpolate each CDF
for i_cdf = 1:n_cdfs
    x = bin_edges{i_cdf};
    v = cdfs{i_cdf};
    cdf_values_interp(:,i_cdf) = interp1(x(1:end-1),v,bin_edges_interp);
    cdf_values_interp(bin_edges_interp < x(1),i_cdf) = 0;
    cdf_values_interp(bin_edges_interp >= x(end),i_cdf) = 1;
end