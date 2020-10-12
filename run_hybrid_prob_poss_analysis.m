% Runscript for hybrid probabilistic-possibilistic uncertainty analysis 
% applied to a socio-economic cost-benefit analysis, accounting for 
% uncertainties due to high-impact low-probability (HILP) events.

%% Set up input parameters etc.

% Number of alpha cuts (including alpha = 0 and alpha = 1)
n_cuts = 4;

% Number of latin hypercube samples
n_lh_samples = 120;

% True if using surrogate modelling option for calculating possibility distributions
use_surrogate_poss_distr = false;

% True if using surrogate modelling option for calculating p-boxes
use_surrogate_p_box = false;

% True if plotting comparison with vertex optimization for possibility
% distributions
do_plot_vertex = false;

% Loading input values for case
input_params_case_study;

% Font sizes
FontSizeLeg = 11;
FontSizeTicks = 11;
FontSizeLabel = 14;

% Set up settings for socio-economic analysis
socioeconomic_analysis_options = set_socioec_options();

% Load parameter values for case study
[parameters_prob, parameters_poss] = input_params_case_study;

% Package settings for the (outer loop of the) hybrid probabilistic-
% possibilistic uncertainty analysis
poss_settings.n_cuts = n_cuts;
poss_settings.n_lh_samples = n_lh_samples;
poss_settings.use_surrogate_poss_distr = use_surrogate_poss_distr;
poss_settings.use_surrogate_p_box = use_surrogate_p_box;

%% Run hybrid probabilistic-possibilistic uncertainty analysis

% Run analysis
[results_poss_distr, results_p_box] = hybrid_prob_poss_analysis(poss_settings,parameters_prob,parameters_poss,socioeconomic_analysis_options);

% Unpacking results for possibility distributions
alpha_vec_poss = results_poss_distr.alpha_vec_poss;
Mean_net_benefits_vec = results_poss_distr.Mean_rel_benefits_vec;
Mean_rel_benefits_vec_vertex = results_poss_distr.Mean_rel_benefits_vec_vertex;
Mean_CB_0_vec = results_poss_distr.Mean_CB_0_vec;
Mean_CB_0_vec_vertex = results_poss_distr.Mean_CB_0_vec_vertex;
Mean_CB_1_vec = results_poss_distr.Mean_CB_1_vec;
Mean_CB_1_vec_vertex = results_poss_distr.Mean_CB_1_vec_vertex;
Prob_positive_net_benefit_vec = results_poss_distr.Prob_positive_net_benefit_vec;

% Unpacking results for p-boxes
bin_edges_Net_benefits = results_p_box.bin_edges_Mean_rel_benefits;
bin_edges_CB = results_p_box.bin_edges_CB;
cdf_values_CB_0 = results_p_box.cdf_values_CB_0;
cdf_values_CB_1 = results_p_box.cdf_values_CB_1;
cdf_values_Net_benefits = results_p_box.cdf_values_Net_benefits;


%% Plot possibility function for the expected net benefits

figure
plot(Mean_net_benefits_vec,alpha_vec_poss,'k','LineWidth',1)
hold on
if do_plot_vertex
    plot(Mean_rel_benefits_vec_vertex,alpha_vec_poss,'k--','LineWidth',1)
end
plot([0 0],[0 1],'k:')
xlim([-1E9 6E9])
xlabel('Expected net benefits of risk-reducing measures (NOK)','Interpreter','latex','FontSize',FontSizeLabel)
ylabel('Possibility distribution $\pi$','Interpreter','latex','FontSize',FontSizeLabel)
set(gca,'TickLabelInterpreter','latex','FontSize',FontSizeTicks)

% Plot uncertainty intervals
fig = figure;
hold on
plot([min(Mean_net_benefits_vec) max(Mean_net_benefits_vec)],[1 1],'LineWidth',1,'Color','k','LineStyle','-');
plot([min(Mean_net_benefits_vec) min(Mean_net_benefits_vec)],[0.9 1.1],'LineWidth',1,'Color','k','LineStyle','-');
plot([max(Mean_net_benefits_vec) max(Mean_net_benefits_vec)],[0.9 1.1],'LineWidth',1,'Color','k','LineStyle','-');
plot([Mean_net_benefits_vec(n_cuts) Mean_net_benefits_vec(n_cuts)],[-1 1],'LineWidth',0.8,'Color','k','LineStyle','none','Marker','x');
plot([0 0],[0 2],'k:')
xlabel('Expected net benefits of risk-reducing measures (NOK)','Interpreter','latex','FontSize',FontSizeLabel)
ylim([0 2])
xlim([-1E9 6E9])
set(gca,'YTick',[])
set(gca,'TickLabelInterpreter','latex','FontSize',FontSizeTicks)
fig.Position(4) = 150;

%% Plot the expected socio-economic costs of alternatives

% Plot possibility distributions
figure
plot(Mean_CB_0_vec,alpha_vec_poss,'r','LineWidth',1)
hold on
if do_plot_vertex
    plot(Mean_CB_0_vec_vertex,alpha_vec_poss,'r--','LineWidth',0.5)
end
plot(Mean_CB_1_vec,alpha_vec_poss,'b','LineWidth',1)
if do_plot_vertex
    plot(Mean_CB_1_vec_vertex,alpha_vec_poss,'b--','LineWidth',0.5)
end

xlim([0 6E9])
xlabel('Expected socio-economic costs of alternatives (NOK)','Interpreter','latex','FontSize',FontSizeLabel)
ylabel('Possibility distribution $\pi$','Interpreter','latex','FontSize',FontSizeLabel)
set(gca,'TickLabelInterpreter','latex','FontSize',FontSizeTicks)
leg = legend({'A: Zero alternative','B: Risk-reducing measures (additional transmission line)'},'EdgeColor','none');
if do_plot_vertex
    leg = legend({'A: Zero alternative ($n_\mathrm{LHS} = 80$)',...
        'A: Zero alternative ($n_\mathrm{LHS} = 0$)',...
        'B: With additional transmission line ($n_\mathrm{LHS} = 80$)',...
        'B: With additional transmission line ($n_\mathrm{LHS} = 0$)',},'EdgeColor','none');
end
set(leg,'Interpreter','LaTeX')
set(leg,'FontSize',FontSizeLeg)

%% Plot probability of having net benefit

figure
plot(Prob_positive_net_benefit_vec,alpha_vec_poss,'k','LineWidth',1)
xlabel('Probability that net benefits $>$ costs for risk-reducing measures','Interpreter','latex','FontSize',FontSizeLabel)
ylabel('Possibility distribution $\pi$','Interpreter','latex','FontSize',FontSizeLabel)
xlim([0 1])
ylim([0 1])

% Plot uncertainty intervals
fig = figure;
hold on
plot([min(Prob_positive_net_benefit_vec) max(Prob_positive_net_benefit_vec)],[1 1],'LineWidth',1,'Color','k','LineStyle','-');
plot([min(Prob_positive_net_benefit_vec) min(Prob_positive_net_benefit_vec)],[0.9 1.1],'LineWidth',1,'Color','k','LineStyle','-');
plot([max(Prob_positive_net_benefit_vec) max(Prob_positive_net_benefit_vec)],[0.9 1.1],'LineWidth',1,'Color','k','LineStyle','-');
plot([Prob_positive_net_benefit_vec(n_cuts) Prob_positive_net_benefit_vec(n_cuts)],[-1 1],'LineWidth',0.8,'Color','k','LineStyle','none','Marker','x');
xlabel('Probability that net benefits > costs for risk-reducing measures','Interpreter','latex','FontSize',FontSizeLabel)
ylim([0 3])
xlim([0 1])
set(gca,'YTick',[])
set(gca,'TickLabelInterpreter','latex','FontSize',FontSizeTicks)
fig.Position(4) = 150;

%% Plot p-boxes

% Plot P-box of net benefit
figure
plot([0 0],[0 1],'k:')

p_net_benefit = plot_p_box(cdf_values_Net_benefits,bin_edges_Net_benefits,[112 48 160]/255);

% Plot best-guess expected value of net benefits
Mean_rel_benefits_0 = Mean_net_benefits_vec(n_cuts);
plot([1 1]*Mean_rel_benefits_0,[0 1],'Color',[112 48 160]/255,'LineWidth',0.7,'LineStyle','--')

set(gca,'TickLabelInterpreter','latex','FontSize',12)
xlabel('Net benefits of risk-reducing measures $\Delta$TC (NOK)','Interpreter','latex','FontSize',FontSizeLabel)
%ylim([0.95 1])  % Zoomed-in view
ylim([0 1])   % Zoomed-out view
xlim([-1E9 6E9])
set(gca,'TickLabelInterpreter','latex','FontSize',FontSizeTicks)


% Plot P-boxes for cost-benefit for each alternative
figure
plot([0 0],[0 1],'k:')

p_CB_0 = plot_p_box(cdf_values_CB_0,bin_edges_CB,'r');
p_CB_1 = plot_p_box(cdf_values_CB_1,bin_edges_CB,'b');

xlabel('Total socio-economic costs TC (NOK)','Interpreter','latex','FontSize',FontSizeLabel)
%ylim([0.95 1])  % Zoomed-in view
ylim([0 1])   % Zoomed-out view
%    plot([0 0],[0 1],'k:')
xlim([-1E9 6E9])
set(gca,'TickLabelInterpreter','latex','FontSize',FontSizeTicks)
leg = legend([p_CB_0,p_CB_1],{'A: Zero alternative','B: Risk-reducing measures (additional transmission line)'});
set(leg,'EdgeColor','none')
set(leg,'Color','none')
set(leg,'Interpreter','LaTeX')
set(leg,'FontSize',FontSizeLeg)

