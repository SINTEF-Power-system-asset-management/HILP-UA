function compare_cdfs(data,alpha,do_plot_log,str_xlabel,str_legend,do_plot_CVaR,cmap)
% COMPARE_CDFS Plots comparison between several cumulative distribution functions
%
% INPUTS:
%   data:
%   Cell array with data for the CDFs, with each element being either a 
%   histogram object or a numerical array
%
%   alpha: 
%   alpha value for Value at Risk / percentile calculation (giving the
%   (1 - alpha) x 100 % percentile)
%
%   do_plot_log: 
%   True if plotting in a log-log plot
%
%   str_xlabel: 
%   Char string for x axis label
%
%   str_legend:
%   Cell array of char strings for the legend
%
%   do_plot_CVaR:
%   Logical being true if plotting Conditional Value-at-Risk (optional;
%   default: true)
%
%   cmap:
%   Colour map (see function COLORMAP; optional)


% Font sizes
FontSizeLeg = 11;
FontSizeTicks = 11;
FontSizeLabel = 13;

% Factor to multiply with colour to make it darker, [0,1], with 1 making it
% unchanged and 0 making all lines black...
factor_darker = 0.75;

% Plotting CVaR by default
if nargin < 6
    do_plot_CVaR = true;
end

% Choose colourmap 
if nargin < 7
    cmap = colormap('cool');
    cmap = cmap .* factor_darker;
end

if ~iscell(data)
    error('Input argument data needs to be a cell array')
end

n_data = numel(data);

mean_alt = zeros(n_data,1);
VaR_alt = zeros(n_data,1);
CVaR_alt = zeros(n_data,1);
h_alt = cell(n_data,1);


for i = 1:n_data
    if isobject(data{i}) && strcmp(class(data{i}),'matlab.graphics.chart.primitive.Histogram')
        h_alt{i} = data{i};
        mean_alt(i) = hist_mean(h_alt{i});
    elseif isnumeric(data{i})
        if i == 1
            fig1 = figure;
            axes1 = axes('Parent',fig1);
        else
            hold on
        end
        h_alt{i} = histogram(data{i});
        mean_alt(i) = mean(data{i});
    else
        error('Unsupported data format for input argument 1')
    end
end

% Calculate (Conditional) Value at Risk
for i = 1:n_data
    [CVaR_alt(i), VaR_alt(i)] = calc_CVaR(data{i},alpha);
end

% Choose one colour for each scenario from colourmap
colours = cell(n_data,1);
if n_data == 1
    colours{i} = cmap(1,:);
else    
    for i = 1:n_data
        i_colour =  round((size(cmap,1)-1) * (i-1)/(n_data-1)) + 1;
        colours{i} = cmap(i_colour,:);
    end
end

%% Plotting
if isobject(data{i}) && strcmp(class(data{i}),'matlab.graphics.chart.primitive.Histogram')
    fig1 = figure;
    axes1 = axes('Parent',fig1);
end

x_alt = cell(n_data,1);
y_alt = cell(n_data,1);
x_max = Inf;
for i = 1:n_data
    [h_alt{i} , h_stairs(i)] = format_cdf(h_alt{i},true,colours{i});
    x_alt{i} = h_stairs(i).XData;
    if x_alt{i}(end) < x_max
        x_max = x_alt{i}(end);
    end
    if do_plot_log
        y_alt{i} = -log10 (1 - h_stairs(i).YData);
    else
        y_alt{i} = h_stairs(i).YData;
    end
end
hold(axes1,'off');

for i = 1:n_data
    h_stairs(i) = plot(axes1,x_alt{i},y_alt{i},'LineWidth',1,'Color',colours{i});
    if i == 1
        hold(axes1,'on');
    end
end
xlabel(str_xlabel,'Interpreter','latex','FontSize',FontSizeLabel)

if do_plot_log    
    x_min = Inf;
    x_max = -Inf;
    y_max = 5;
    i_start_alt = zeros(n_data,1);
    for i = 1:n_data
        x_min = min(0.1 * mean_alt(i),x_min);
    end
    for i = 1:n_data
        x_max = max(x_max, x_alt{i}(end));
        i_start_alt(i) = find(x_alt{i} > 0,1);
        plot(axes1,[x_min x_alt{i}(i_start_alt(i))],[y_alt{i}(i_start_alt(i)) y_alt{i}(i_start_alt(i))],'LineWidth',1,'Color',colours{i});       
    end
    set(axes1,'XScale','log')
    xlim([x_min, x_max])
    set(axes1,'YTick',[0 1 2 3 4],'YTickLabel',{'0','0.9','0.99','0.999','0.9999'});
    
    ylim([0 y_max])
else
    y_max = 1;
    xlim([0 x_max])
    %ylim([0.95 y_max])
    ylim([0 y_max])
end

ylabel('Cumulative distribution function $F$','Interpreter','latex','FontSize',FontSizeLabel)
set(axes1,'TickLabelInterpreter','latex','FontSize',FontSizeTicks)
set(axes1,'Box','off')

for i = 1:n_data
    plot([mean_alt(i) mean_alt(i)],[0 y_max],'--','Color',colours{i})    
    plot([VaR_alt(i) VaR_alt(i)],[0 y_max],':','Color',colours{i})
    if do_plot_CVaR
        plot([CVaR_alt(i) CVaR_alt(i)],[0 y_max],'-.','Color',colours{i})
    end
end

leg = legend(h_stairs,str_legend);
set(leg,'EdgeColor','none','Location','best')
set(leg,'Interpreter','LaTeX')
set(leg,'FontSize',FontSizeLeg)