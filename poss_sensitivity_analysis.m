function x_break_even = poss_sensitivity_analysis(inputs_poss_sensitivity)
% POSS_SENSITIVITY_ANALYSIS Runs possibilistic sensitivity analysis
%   
%   X_BREAK_EVEN = POSS_SENSITIVITY_ANALYSIS(INPUTS_POSS_SENSITIVITY) runs 
%   a sensitivity analysis considering the sensitivities in the conclusions of a
%   cost-benefit analysis considering changes in parameters with epistemic
%   uncertainties (represented possibilistically). Currently it only does a
%   mapping of the epistemic uncertainty space and calculates the break-even
%   values for the uncertain parameters X_BREAK_EVEN for which the costs of 
%   risk-reducing measures balance the benefits. (TODO: Implement calculation 
%   and plotting of sensitivity indices.)
%
% See also: HYBRID_PROB_POSS_ANALYSIS.M

% TODO: Put plotting functionality in a separate function or runscript?

%% Parameters for sensitivity analysis

% Number of samples per dimension when sampling the response surface
n_samples_per_dim = 50;

% Number of samples in search direction (along each dimension) when
% searching for the break-even value of the cost-benefit analysis
n_samples_search = 1000;

% True if only plotting one (two-dimensional) slice of the multi-dimensional 
% epistemic uncertainty space
only_plot_2D = true;

% Labels for uncertain parameters (1 - Annual frequency of HILP events,
% 2 - Residual risk factor, 3 - Annual load growth, 4 - Annual SoS value growth)
param_labels = {'$\lambda$ (1/year)','$f_\mathrm{red}$','$f_P$','$f_c$'};
param_descr = {'Annual frequency of HILP events','Residual risk factor','Annual load growth','Annual SoS value growth'};

% Which dimensions to plot if only plotting one slice
combs_param_2D = [1 4];

% True if only plotting samples corresponding to vertices of
% hyperrectangles corresponding to alpha-cuts; false if also plotting
% samples from Latin Hypercube Sampling
do_plot_only_vertex_samples = false;

% Text strings for legend in mapping of uncertainty space
str_legend = {'Negative expected net benefit','Positive expected net benefit','For best-guess parameter values $x_0$'};

% Font sizes
FontSizeLeg = 11;
FontSizeTicks = 11;
FontSizeLabel = 14;

% Colours for plotting
color_negative = [220,180,180]/255;
color_positive = [70,100,70]/255;
color_best_guess = [255,0,0]/255;
color_negative_light = [250,240,240]/255;
color_positive_light = [190,210,190]/255;

%% Pre-processing

% Unpacking input parameters from hybrid probabilistic-possibilistic
% uncertainty analysis
Mean_net_benefits_samples = inputs_poss_sensitivity.Mean_net_benefits_samples;
x_poss_sample = inputs_poss_sensitivity.x_poss_sample;
idx_samples_vertex = inputs_poss_sensitivity.idx_samples_vertex;
x_best_guess = inputs_poss_sensitivity.x_best_guess;

% Number of samples from the epistemic uncertainty space
n_samples = size(x_poss_sample,1);

% Number of parameters with uncertainty represented possibilistically
n_params_poss = size(x_poss_sample,2);

% Upper and lower bounds of uncertainty hyperrectangle
x_lower = min(x_poss_sample);
x_upper = max(x_poss_sample);

%% Plot before applying surrogate modelling

% This is only used for comparison with results after applying surrogate
% modelling if the uncertainty space mapping is illustrated for a single
% 2D slice of the uncertainty space
if only_plot_2D
    
    % Samples to plot
    if do_plot_only_vertex_samples
        idx_samples_sel = idx_samples_vertex;
    else
        idx_samples_sel = true(n_samples,1);
    end
    
    % Find those samples that are "close" to the best-guess value for the
    % dimensions of the uncertain parameter space that are not included in
    % the two-dimensional slice
    idx_close_to_best_guess = true(n_samples,1);
    for i_param = 1:n_params_poss
        if ~ismember(i_param,combs_param_2D)
            max_dist_close = 0.2 * ( x_upper(1,i_param) - x_lower(1,i_param) );
            idx_close_to_best_guess = idx_close_to_best_guess & abs(x_poss_sample(:,i_param) - x_best_guess(1,i_param)) < max_dist_close;
        end
    end

    % Uncertain parameters to plot slice for
    i_param_x = combs_param_2D(1);
    i_param_y = combs_param_2D(2);
    
    % Find sample points with positive expected net benefit
    idx_samples_positive = (Mean_net_benefits_samples > 0);
    
    % Plot mapping of epistemic uncertainty space (without using surrogate
    % modelling for the expected net benefit); the samples that are far
    % away from the 2D slice in the multi-dimensional parameter space are
    % plotted with lighter colours
    figure
    title('Sensitivity analysis: Mapping of uncertainty space (without surrogate modelling)')
    hold on
    plot(x_poss_sample(~idx_samples_positive & idx_samples_sel & idx_close_to_best_guess,i_param_x),x_poss_sample(~idx_samples_positive & idx_samples_sel & idx_close_to_best_guess,i_param_y),'LineStyle','none','Color',color_negative,'Marker','.','MarkerSize',15)
    plot(x_poss_sample(idx_samples_positive & idx_samples_sel & idx_close_to_best_guess,i_param_x),x_poss_sample(idx_samples_positive & idx_samples_sel & idx_close_to_best_guess,i_param_y),'LineStyle','none','Color',color_positive,'Marker','.','MarkerSize',15)
    plot(x_best_guess(i_param_x),x_best_guess(i_param_y),'LineStyle','none','Color',color_best_guess,'Marker','x','MarkerSize',10,'LineWidth',1)
    plot(x_poss_sample(~idx_samples_positive & idx_samples_sel & ~idx_close_to_best_guess,i_param_x),x_poss_sample(~idx_samples_positive & idx_samples_sel & ~idx_close_to_best_guess,i_param_y),'LineStyle','none','Color',color_negative_light,'Marker','.','MarkerSize',15)
    plot(x_poss_sample(idx_samples_positive & idx_samples_sel & ~idx_close_to_best_guess,i_param_x),x_poss_sample(idx_samples_positive & idx_samples_sel & ~idx_close_to_best_guess,i_param_y),'LineStyle','none','Color',color_positive_light,'Marker','.','MarkerSize',15)
    xlim([ 0.9*min(x_poss_sample(:,i_param_x)) 1.1*max(x_poss_sample(:,i_param_x)) ])
    ylim([ 0.9*min(x_poss_sample(:,i_param_y)) 1.1*max(x_poss_sample(:,i_param_y)) ])
    xlabel(param_labels{i_param_x},'Interpreter','latex','FontSize',FontSizeLabel)
    ylabel(param_labels{i_param_y},'Interpreter','latex','FontSize',FontSizeLabel)
    leg = legend(str_legend,'Location','Best');
    set(leg,'Interpreter','LaTeX')
    set(leg,'FontSize',FontSizeLeg)
end

%% Make surrogate model for the expected net benefit and plot

% Design matrix for quadratic multiple regression 
% (NB: needs statistics toolbox)
x_design = x2fx(x_poss_sample,'quadratic');

% Finding regression coefficients for quadratic response surface surrogate
% model
beta = (x_design' * x_design)\ (x_design' * Mean_net_benefits_samples);

% We will combine two and two of the uncertain parameters to make a series
% of 2D plots ("slices" through the uncertainty hyperrectangle)
if ~only_plot_2D
    combs_param = nchoosek(1:n_params_poss,2);
    n_params_plot = n_params_poss;
else
    combs_param = combs_param_2D;
    n_params_plot = 2;
end
n_combs_param = size(combs_param,1);

% Plot mapping of epistemic uncertainty space (without using surrogate
% modelling for the expected net benefit)
figure
for i_comb = 1:n_combs_param
    i_param_x = combs_param(i_comb,1);
    i_param_y = combs_param(i_comb,2);
    if n_params_plot > 2
        pos = (i_param_y-2) * (n_params_poss-1) + i_param_x;
        subplot(n_params_poss-1,n_params_poss-1,pos)
    end
        
    % We will fix all other parameters than those combined to their
    % best-guess values
    idx_param_fixed = ~ismember(1:n_params_poss,combs_param(i_comb,:));
    x_fixed = NaN(1,n_params_poss);
    x_fixed(idx_param_fixed) = x_best_guess(idx_param_fixed);
    
    % Sample points on the response surface using a linearly spaced grid
    x_fit = select_grid_samples(n_samples_per_dim,x_lower,x_upper,x_fixed);    
    x_fit_design = x2fx(x_fit,'quadratic');
    Mean_net_benefits_fit = x_fit_design * beta;
    
    % The sampled points with positive net benefit
    idx_positive_fit = (Mean_net_benefits_fit > 0);
    
    % Plotting
    plot(x_fit(~idx_positive_fit,i_param_x),x_fit(~idx_positive_fit,i_param_y),'LineStyle','none','Color',color_negative,'Marker','.','MarkerSize',15)
    hold on
    plot(x_fit(idx_positive_fit,i_param_x),x_fit(idx_positive_fit,i_param_y),'LineStyle','none','Color',color_positive,'Marker','.','MarkerSize',15)
    plot(x_best_guess(i_param_x),x_best_guess(i_param_y),'LineStyle','none','Color',color_best_guess,'Marker','x','MarkerSize',10,'LineWidth',1)
    
    xlim([ 0.9*min(x_fit(:,i_param_x)) 1.1*max(x_fit(:,i_param_x)) ])
    ylim([ 0.9*min(x_fit(:,i_param_y)) 1.1*max(x_fit(:,i_param_y)) ])
    xlabel(param_labels{i_param_x},'Interpreter','latex','FontSize',FontSizeLabel)
    ylabel(param_labels{i_param_y},'Interpreter','latex','FontSize',FontSizeLabel)
    set(gca,'TickLabelInterpreter','latex','FontSize',FontSizeTicks)
    set(gca,'Box','off')
    
    if i_comb == 1
        leg = legend(str_legend,'Location','Best');
        set(leg,'Interpreter','LaTeX')
        set(leg,'FontSize',FontSizeLeg)
        set(leg,'EdgeColor','none')
        ylims = get(gca,'YLim');
        ylim([ylims(1) 1.1 * ylims(2)])
    end
end

if only_plot_2D
    title('Sensitivity analysis: Mapping of uncertainty space (with surrogate modelling)')
end

%% Searching for the break-even value of the cost-benefit analysis

% Initialize parameter vector for searching for break-even
x_search_init = repmat(x_best_guess,n_samples_search,1);

% Initialize parameter values for break-even
x_break_even = NaN(1,n_params_poss);

disp('Break-even values of uncertain parameters in the cost-benefit analysis:')

% Search for the parameter value of each uncertain parameter where the
% cost-benefit analysis breaks even, assuming the other parameters to be
% fixed to the best-guess value
for i_param = 1:n_params_poss
    x_search = x_search_init;
    
    % Vector with all parameters except the search direction set to best-guess value
    x_search_i = linspace(x_lower(i_param),x_upper(i_param),n_samples_search)';
    x_search(:,i_param) = x_search_i;
    
    % Calculate mean net benefit along search direction
    x_search_design = x2fx(x_search,'quadratic');
    Mean_net_benefits_search = x_search_design * beta;
    
    % Finding break-even point along search direction where net benefits become positive
    x_break_even_i = x_search_i(find(Mean_net_benefits_search > 0,1,'first'));
    if ~isempty(x_break_even_i)
        x_break_even(i_param) = x_search_i(find(Mean_net_benefits_search > 0,1,'first'));
    end
    
    disp([param_descr{i_param} ': ' num2str(x_break_even(i_param))])
end
