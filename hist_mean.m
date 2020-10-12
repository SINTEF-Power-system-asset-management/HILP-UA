function y_mean = hist_mean(h)
% HIST_MEAN Calculates mean value from a histogram
%
%   Y_MEAN = HIST_MEAN(H) returns the (approximate) mean value Y_MEAN for 
%   the data that the histogram object H is based on.

if ~isobject(h) || ~strcmp(class(h),'matlab.graphics.chart.primitive.Histogram')
    error('Input argument needs to be a histogram')
end
   
h.Normalization = 'pdf';
BinCenters = ( h.BinEdges(1:end-1) + h.BinEdges(2:end) )/2;
y_mean = sum(BinCenters .* h.BinWidth .* h.Values);
