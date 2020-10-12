function IC = interr_cost_model(l,P,r,c,c_r)
% INTERR_COST_MODEL Interruption cost model
% 
%   IC = INTERR_COST_MODEL(L,P,R,C) returns the power supply
%   interruption costs IC for frequency of interruption events L (\lambda,
%   events per year), interrupted power P (MW), interruption duration R 
%   (hours), specific interruption cost C (NOK/MWh). The interruption costs
%   are given in units of NOK/year; however, if L = 1, the results can be
%   interpreted as the interruption costs of a single HILP event (given
%   that it occurs).
%
%   IC = INTERR_COST_MODEL(L,P,R,C,C_R) returns the interruption costs
%   accounting for a hypothetical linear relationship between the
%   interruption duration and the specific interruption costs. (This is
%   used as a simple test case for correlations between input parameters.)

if nargin < 5
    c_r = 0;
end

IC = l .* P .* r .* c .* (1 + c_r*r);