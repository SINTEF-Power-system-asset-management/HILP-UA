function h_p_box = plot_p_box(cdf_values, bin_edges,colour)
% PLOT_P_BOX Plots a P-box for of aleatory and epistemic uncertainty
%
%   PLOT_P_BOX(CDF_VALUES, BIN_EDGES,I_CDF_CRISP,COLOUR) plots the probability 
%   bounds or the "P-box" of CDFs illustrating the contributions of aleatory and 
%   epistemic uncertainty defined by a set of N_CDF = 2 or N_CDF = 3 
%   Cumulative Distribution Functions (CDFs, defined for a set of histogram 
%   bin edges BIN_EDGES). CDFS is a N_BINS x N_CDF array and N_BINS is a 
%   N_BINS x 1 or N_BINS x N_CDF array. The first column in CDF_VALUES is
%   the lower probability bound and the second column is the upper
%   probability bound; the third column (if present) is the "crisp" CDF
%   corresponding to the best-guess (most likely) probability distribution.

if nargin < 3
    colour = [0.25 0.25 0.25];
end

% True if using logarithmic x axis
use_log_x = false;

hold on

cdf_lower = cdf_values(:,1);
cdf_upper = cdf_values(:,2);

if size(cdf_values,2) == 2
    do_plot_crisp = false;
elseif size(cdf_values,2) == 3
    do_plot_crisp = true;
    cdf_crisp = cdf_values(:,3);
else
    error('Only 2 or 3 columns supported for input argument cdf_values')
end

% Set up for plotting bounds
x = repmat(bin_edges',2,1);
x = x(2:end);

y_upper = repmat(cdf_upper',2,1);
y_upper = y_upper(1:end-1);

y_lower = repmat(cdf_lower',2,1);
y_lower = y_lower(1:end-1);

% Plotting
h_p_box = fill([x fliplr(x)],[y_upper fliplr(y_lower)],colour);
set(h_p_box,'LineWidth',0.5,'FaceAlpha',0.5,'EdgeColor',colour)
hold on
if do_plot_crisp
    y_crisp = repmat(cdf_crisp',2,1);
    y_crisp = y_crisp(1:end-1);
    h_crisp = plot(x,y_crisp,'LineWidth',1,'Color',colour,'LineStyle','-');
end

ylabel('Cumulative distribution function $F$','Interpreter','latex','FontSize',13)

set(gca,'Box','on')

if use_log_x
    set(gca,'XScale','log');
end

if do_plot_crisp
    leg = legend([h_p_box,h_crisp],{'Probability bounds (between upper and lower CDFs)','"Best guess" cumulative distribution function'});
else
    leg = legend([h_p_box,h_crisp],{'Probability bounds (between upper and lower CDFs)'});
end
set(leg,'EdgeColor','none')
set(leg,'Color','none')
set(leg,'Interpreter','LaTeX')
set(leg,'FontSize',11)