# HILP-UA: HILP Uncertainty Analysis

MATLAB code for hybrid probabilistic-possibilistic uncertainty analysis 
applied to socio-economic cost-benefit analysis for power system development
decisions involving high-impact low-probability (HILP) events.

Source code used to generate the results reported in the manuscript:
I. B. Sperstad, G. H. Kj�lle, and E. �. Norum, "Accounting for uncertainties 
due to high-impact low-probability events in power system development", 
(Submitted for peer review), 2020.

The code base includes test case data that can be used to generate the same 
type of results as reported in the manuscript. However, since the case 
study considered a real case in Norway and the load time series that was
used is confidential, a synthetic load time series is included that
produces similar results as in the manuscript but do not reproduce the
results exactly.

The code requires the MATLAB Optimization Toolbox and the MATLAB 
Statistics Toolbox to run if the optional surrogate model functionality is 
enabled.

The script run_prob_analysis_single_event.m runs a probabilisic uncertainty
analysis for a single HILP event (given that it occurs) and can be run to
to (approximately) reproduce Figure 3 in the manuscript.

The script run_prob_socioec_analysis.m runs a probabilistic socio-economic 
cost-benefit analysis considering interruption costs due to HILP events and 
can be run to (approximately) reproduce Figures 4 and 5 in the manuscript.

The script run_hybrid_prob_poss_analysis.m runs a hybrid probabilistic-
possibilistic uncertainty analysis applied to the socio-economic cost-benefit 
analysis, accounting for uncertainties due to HILP events and can be used to 
reproduce Figures 6, 7 and 8 in the manuscript.

The input data of the case study are specified in the function 
input_params_case_study.m. The options for the socio-economic cost-benefit 
analysis can be set in the function set_socioec_options.m.

Author: Iver Bakken Sperstad, SINTEF Energy Research 
(iver.bakken.sperstad@sintef.no)