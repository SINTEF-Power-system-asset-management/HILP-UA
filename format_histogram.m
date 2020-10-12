function format_histogram(h,incl_label)
% FORMAT_CDF Formats plot of a histogram
% 
%   FORMAT_CDF(H,INCL_LABEL) formats an input histogram object handle H. 
%   An y-axis label is printed if INCL_LABEL = true.

set(h,'Normalization','pdf')
set(h,'EdgeColor','none')
set(h,'FaceColor',[0.5 0.5 0.5])
hold on
stairs([h.BinEdges(1)-h.BinWidth h.BinEdges],[0 h.Values 0],'k')

if nargin > 1 && incl_label == true
    ylabel('Estimated probability density function')
end