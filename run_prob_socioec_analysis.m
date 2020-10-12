% Runscript for a probabilistic socio-economic cost-benefit analysis
% considering interruption costs due to high-impact low-probability (HILP) 
% events

%% Settings

% True if considering a portfolio perspective with n_projects identical
% projects
do_consider_portfolio = false;
n_projects = 20;

% True if plotting only zero alternative in interruption cost plots with
% linear axes
do_plot_only_zero = false;

% True if plotting results (overriding settings in set_socioec_options.m)
do_plot = true;

% Font sizes
FontSizeLeg = 11;
FontSizeTicks = 11;
FontSizeLabel = 14;

%% Load input parameter values

% Load input parameter values for case study
[parameters_prob, parameters_poss] = input_params_case_study;

% Percentage increase in valuation of lost load per year (best-guess value)
r_incr_c = parameters_poss.r_incr_c_input.core;

% Expected percentage-wise load growth per year (best-guess value)
r_load_growth = parameters_poss.r_load_growth_input.core;

% Fraction of irreducible HILP events (best-guess value)
f_irred = parameters_poss.f_irred_input.core;

% Frequency of occurrence of HILP events (best-guess value)
l = parameters_poss.l_input.core;

%% Packaging input parameters for socio-economic analysis

% Set options to use for socio-economic analysis
socioeconomic_analysis_options = set_socioec_options();
socioeconomic_analysis_options.do_plot = do_plot;
do_show_BC_and_not_CB = socioeconomic_analysis_options.do_show_BC_and_not_CB;
n_MC_sims = socioeconomic_analysis_options.n_MC_sims;

% Parameters treated as having fixed values in this function (but that may have epistemic uncertainty)
parameters_fixed.r_load_growth = r_load_growth;
parameters_fixed.r_incr_c = r_incr_c;
parameters_fixed.f_irred = f_irred;
parameters_fixed.l = l;

if do_consider_portfolio
    parameters_prob.C_0_input.mode = parameters_prob.C_0_input.mode * n_projects;
    parameters_prob.C_0_input.lower = parameters_prob.C_0_input.lower * n_projects;
    parameters_prob.C_0_input.upper = parameters_prob.C_0_input.upper * n_projects;
    parameters_prob.C_1_input.mode = parameters_prob.C_1_input.mode * n_projects;
    parameters_prob.C_1_input.lower = parameters_prob.C_1_input.lower * n_projects;
    parameters_prob.C_1_input.upper = parameters_prob.C_1_input.upper * n_projects;
    parameters_fixed.l = parameters_fixed.l * n_projects;
    socioeconomic_analysis_options.do_allow_max_1_event = false;
end

%% Run probabilistic socio-economic cost-benefit analysis

[IC_alt_0,IC_alt_1,CB_0,CB_1,Rel_benefits] = prob_socioec_analysis(parameters_prob,parameters_fixed,socioeconomic_analysis_options);

IC_alt_0_mean = mean(IC_alt_0);
ste_IC_alt_0 = std(IC_alt_0)/sqrt(n_MC_sims);

IC_alt_1_mean = mean(IC_alt_1);
ste_IC_alt_1 = std(IC_alt_1)/sqrt(n_MC_sims);

CB_alt_0_mean = mean(CB_0);
ste_CB_alt_0 = std(CB_0)/sqrt(n_MC_sims);

CB_alt_1_mean = mean(CB_1);
ste_CB_alt_1 = std(CB_1)/sqrt(n_MC_sims);

Mean_rel_benefits = mean(Rel_benefits);
steRel_benefits = std(Rel_benefits)/sqrt(n_MC_sims);


%% Plotting interruption cost results in zoomed out mode

alpha = 0.001;
do_plot_log = false;
do_plot_CVaR = false;
str_xlabel = 'Interruption costs for grid development alternative (NOK)';
str_legend = {'A: Zero alternative','B: Risk-reducing measures (additional transmission line)'};
cmap = [1 0 0; 0 0 1];
if do_plot_only_zero
    compare_cdfs({IC_alt_0},alpha,do_plot_log,str_xlabel,str_legend{1},do_plot_CVaR,cmap)
else
    compare_cdfs({IC_alt_0,IC_alt_1},alpha,do_plot_log,str_xlabel,str_legend,do_plot_CVaR,cmap)
end
set(gca,'Box','on')
xlim([0 5]*1E9)
ylim([0 1])

% Printing estimated VaR and CVaR 
[CVaR_IC_alt_0, VaR_IC_alt_0] = calc_CVaR(IC_alt_0,alpha)
[CVaR_IC_alt_1, VaR_IC_alt_1] = calc_CVaR(IC_alt_1,alpha)

%% Plotting interruption cost results in zoomed out mode

alpha = 0.001;
do_plot_log = false;
do_plot_CVaR = false;
str_xlabel = 'Interruption costs for grid development alternative (NOK)';
str_legend = {'A: Zero alternative','B: Risk-reducing measures (additional transmission line)'};
cmap = [1 0 0; 0 0 1];
if do_plot_only_zero
    compare_cdfs({IC_alt_0},alpha,do_plot_log,str_xlabel,str_legend{1},do_plot_CVaR,cmap)
else
    compare_cdfs({IC_alt_0,IC_alt_1},alpha,do_plot_log,str_xlabel,str_legend,do_plot_CVaR,cmap)
end
set(gca,'Box','on')
xlim([0 10]*1E9)
ylim([0.95 1])



%% Plotting interruption cost results in log-log plot for both grid development alternatives

alpha = 0.001;
do_plot_log = true;
do_plot_CVaR = false;
str_xlabel = 'Interruption costs for grid development alternative (NOK)';
str_legend = {'A: Zero alternative','B: Risk-reducing measures (additional transmission line)'};
cmap = [1 0 0; 0 0 1];
compare_cdfs({IC_alt_0,IC_alt_1},alpha,do_plot_log,str_xlabel,str_legend,do_plot_CVaR,cmap)
set(gca,'Box','on')


%% Plotting net benefits of risk-reducing measures

if do_plot
    fig1 = figure;
    axes1 = axes('Parent',fig1);
    h_benefits = histogram(Rel_benefits);
    
    % In case of inaccurate binning for the lowest values, reset the
    % first bin edge so that it is not lower than the lowest result value
    bin_edges = h_benefits.BinEdges;
    bin_edges(1) = min(Rel_benefits) - eps;
    set(h_benefits,'BinEdges',bin_edges)
    
    [h_benefits, h_stairs_benefits] = format_cdf(h_benefits,true,[112 48 160]/255);
    hold(axes1,'off');
    
    x_benefits = h_stairs_benefits.XData;        
    if do_plot_log
        y_benefits = -log10 (1 - h_stairs_benefits.YData);     
    else
        y_benefits = h_stairs_benefits.YData;
    end
    
    h_stairs_benefits = plot(axes1,x_benefits,y_benefits,'LineWidth',1,'Color',[112 48 160]/255);
    hold(axes1,'on');
    h_mean_benefit = plot(axes1,[Mean_rel_benefits Mean_rel_benefits],[0 1],'--','Color',[112 48 160]/255);   
    plot([0 0],[0 1],'k:')   
    xlabel('Net benefits of risk-reducing measures (NOK)')
    ylabel('Cumulative distribution function')
    
    if do_plot_log
        y_max = 5;
        set(axes1,'XScale','log')
        xlim([x_benefits(find(x_benefits>0,1)) x_benefits(end)])
        set(axes1,'YTick',[0 1 2 3 4],'YTickLabel',{'0','0.9','0.99','0.999','0.9999'});
        ylim([0 y_max])        
    else
        y_max = 1;   
        xlim([-1E9 5E9])
        ylim([0.95 y_max])
    end
    
    leg = legend([h_stairs_benefits,h_mean_benefit],{'Cumulative distribution function','Expected value'});
    set(leg,'EdgeColor','none')
    set(leg,'Color','none')
    set(leg,'Location','best')
    set(leg,'Interpreter','LaTeX')
    set(leg,'FontSize',10)

end

%% Plot total socio-economic costs of grid development alternatives

if do_plot
    figure
    h_CB_alt_0 = histogram(CB_0);
    
    % In case of inaccurate binning for the lowest values, reset the
    % first bin edge so that it is not lower than the lowest result value
    bin_edges = h_CB_alt_0.BinEdges;
    bin_edges(1) = min(CB_0) - eps;
    set(h_CB_alt_0,'BinEdges',bin_edges)
    
    [h_CB_alt_0, h_stairs_0] = format_cdf(h_CB_alt_0,true,'r');
    
    hold on
    h_CB_alt_1 = histogram(CB_1);
    
    % In case of inaccurate binning for the lowest values, reset the
    % first bin edge so that it is not lower than the lowest result value
    bin_edges = h_CB_alt_1.BinEdges;
    bin_edges(1) = min(CB_1) - eps;
    set(h_CB_alt_1,'BinEdges',bin_edges)
    
    [h_CB_alt_1, h_stairs_1] = format_cdf(h_CB_alt_1,true,'b');       
    if do_show_BC_and_not_CB
        xlabel('Socio-economic impacts of grid development alternative (NOK)')
    else
        xlabel('Total socio-economic costs (NOK)')
    end
    ylabel('Cumulative distribution function')     
    
    plot([0 0],[0 1],'k:')
    plot([CB_alt_0_mean CB_alt_0_mean],[0 1],'--r')
    plot([CB_alt_1_mean CB_alt_1_mean],[0 1],'--b')
    
    if do_show_BC_and_not_CB        
        ylim([0.0 0.05])
        xlim([-IC_alt_0_mean*100 0])
    else
        ylim([0.95 1])
        xlim([-1E9 5E9])
    end
    
    leg = legend([h_stairs_0, h_stairs_1],{'Zero alternative','With risk-reducing measure'});
    set(leg,'EdgeColor','none','Location','best')
end

%% Plot just the expected socio-economic costs of alternatives

if do_plot   
    figure
    hold on
    plot([CB_alt_0_mean CB_alt_0_mean],[0 1],'--r')
    plot([CB_alt_1_mean CB_alt_1_mean],[0 1],'--b')
    
    xlim([0 6E9])
    xlabel('Expected socio-economic costs of alternatives (NOK)','Interpreter','latex','FontSize',FontSizeLabel)
    set(gca,'TickLabelInterpreter','latex','FontSize',FontSizeTicks)
    leg = legend({'A: Zero alternative','B: Risk-reducing measures (additional transmission line)'},'EdgeColor','none');
    set(leg,'Interpreter','LaTeX')
    set(leg,'FontSize',FontSizeLeg)
    set(gca,'YTick',[])
end