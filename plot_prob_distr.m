function plot_prob_distr(IC,do_plot_std,do_plot_VaR,alpha,do_plot_CVaR)
% PLOT_PROB_DISTR Plots probability density function for interruption costs
%
% INPUTS: 
%   do_plot_std: 
%   True if plotting plus and minus one standard deviation
%   
%   do_plot_VaR:
%   Trur if plotting Value-at-Risk at level alpha
%
%   alpha:
%   Level defining Value-at-Risk [0,1]
%
%   do_plot_CVaR:
%   True if plotting Conditional Value-at-Risk

% TODO: Add functionality for colouring high-impact tail to illustrate
% Value-at-Risk


if nargin < 2
    do_plot_std = false;
end

if nargin == 3
    error('If plotting Value-at-Risk, the level parameter \alpha must be provided as the 4th argument')
end

if nargin < 3
    do_plot_VaR = false;
else
    if do_plot_VaR == true && do_plot_std == true
        error('Not supporting plotting both standard deviation and Value-at-Risk')
    end
end

if nargin < 5
    do_plot_CVaR = false;
else
    if do_plot_CVaR && ~ do_plot_VaR
        error('Not supporting plotting Conditional Value-at-Risk without Value-at-Risk')
    end
end

% Defining formatting
FontSizeLeg = 12;
FontSizeTicks = 11;
FontSizeLabel = 14;
linewidth = 1;
linestyle_mean = '--';
linestyle_conf = ':';
linestyle_CVaR = '-.';
linestyle_VaR = ':';
color_lines = 'k';

% Calculating expected interruption cost (NOK)
IC_mean =  mean(IC);
IC_max = max(IC);
IC_std = sqrt(var(IC));

% Plot probability density function
figure
h = histogram(IC,100);
[h, ~] = format_pdf(h,true,'k');

xlabel('Interruption costs due to HILP event (NOK)','FontSize',FontSizeLabel,'Interpreter','latex')

p_IC = h.Values;
IC_bin_start = h.BinEdges;

% Calculate Conditional Value-at-Risk and Value-at-Risk
[CVaR, VaR] = calc_CVaR(IC,alpha);

% Plottng mean value
[~, i_bin_closest] = min(abs(IC_bin_start - IC_mean));
p_IC_mean = p_IC(i_bin_closest);
p_mean = plot([IC_mean IC_mean],[0 p_IC_mean],'Color',color_lines,'LineWidth',linewidth,'LineStyle',linestyle_mean);

% Plot plus and minus one standard deviation
if do_plot_std
    IC_minus_std = IC_mean - IC_std;
    [~, i_bin_closest] = min(abs(IC_bin_start - IC_minus_std));
    p_IC_minus_std = p_IC(i_bin_closest);
    p_minus = plot([IC_minus_std IC_minus_std],[0 p_IC_minus_std],'Color',color_lines,'LineWidth',linewidth,'LineStyle',linestyle_conf);
    
    IC_plus_std = IC_mean + IC_std;
    [~, i_bin_closest] = min(abs(IC_bin_start - IC_plus_std));
    p_IC_plus_std = p_IC(i_bin_closest);    
    p_plus = plot([IC_plus_std IC_plus_std],[0 p_IC_plus_std],'Color',color_lines,'LineWidth',linewidth,'LineStyle',linestyle_conf);
    
    leg = legend([p_mean, p_minus],{'Expected value',' plus/minus one standard deviation'});
end

% Plot Value-at-Risk
if do_plot_VaR
    p_VaR = plot([VaR VaR],[0 p_IC_mean],'LineStyle',linestyle_VaR,'Color',color_lines,'LineWidth',linewidth);   
    
    % Plot Conditional Value-at-Risk
    if do_plot_CVaR
        p_CVaR = plot([CVaR CVaR],[0 p_IC_mean],'LineStyle',linestyle_CVaR,'Color',color_lines,'LineWidth',linewidth)
        leg = legend([p_mean, p_VaR, p_CVaR],{'Expected value','Value-at-Risk','Conditional Variance-at-Risk'});
    else
        leg = legend([p_mean, p_VaR],{'Expected value','Value-at-Risk'});
    end
end

if ~do_plot_std && ~do_plot_VaR
    leg = legend([p_mean],{'Expected value'});
end

set(leg,'EdgeColor','none','FontSize',FontSizeLeg,'Interpreter','latex')
xlim([0 IC_max])
set(gca,'TickLabelInterpreter','latex','FontSize',FontSizeTicks)

