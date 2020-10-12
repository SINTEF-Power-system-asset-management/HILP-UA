function [CVaR, VaR] = calc_CVaR(y,alpha)
% CALC_CVAR Calculate Conditional Value at Risk
%
%   [CVaR, VaR] = CALC_CVAR(y,alpha) returns the Conditional Value at Risk
%   CVaR and optionally the Value at Risk VaR at level alpha for the data 
%   set y (which can be either a numerical array or a histogram object)

    if isobject(y) && strcmp(class(y),'matlab.graphics.chart.primitive.Histogram')
        y.Normalization = 'pdf';
        bin_edges = y.BinEdges;
        pdf_y = y.Values;
        y.Normalization = 'cdf';
        cdf_y = y.Values;
    elseif isnumeric(y)
        [pdf_y,bin_edges] = histcounts(y,'Normalization','pdf');
        [cdf_y] = histcounts(y,'Normalization','cdf');
    else
        error('Unsupported data format for input argument 1')
    end
          
    [~,i_bin_alpha] = min(abs(cdf_y - (1-alpha)));
    VaR = bin_edges(i_bin_alpha);
    bin_width = bin_edges(2)-bin_edges(1);

    % Calculating Conditional Value at Risk
    y_at_risk = zeros(size(pdf_y));
    y_at_risk(i_bin_alpha+1:end) = bin_edges(i_bin_alpha+1:end-1) + bin_width/2;
    CVaR = 1/alpha*(bin_width*pdf_y)*y_at_risk';
