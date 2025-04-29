function plotFittingProgress(res, varargin)
%
% plotFittingProgress(res, varargin)
%
% Description: Plot intermediate values stored throughout DLAG model
%              fitting vs iteration number. Visualizing these values can
%              help with troubleshooting.
%
% Arguments:
%
%     Required:
%
%     res    -- structure containing DLAG fitting results. Contains
%               (relevant) fields
%
%               xDim_across -- int; number of across-group dimensions
%               xDim_within -- (1 x numGroups) array; number of within-
%                              group dimensions for each group  
%               estParams   -- structure containing estimated DLAG model
%                              parameters
%               binWidth    -- float; bin width or sample period, in units 
%                              of time (e.g., ms)
%               iterTime    -- (1 x numIters) array; iterTime(i) contains
%                              the amount of clock time, in sec, EM 
%                              iteration i took to complete.
%               D           -- (1 x numIters) cell array; the estimated 
%                              DLAG delay matrix after each EM iteration.
%               gams_across -- (1 x numIters) cell arry; estimated 
%                              gamma_across after each EM iteration.
%               gams_within -- (1 x numGroups) cell arry;
%                              gams_within(i) -- (1 x numIters) cell array;
%                              estimated gamma_within for group i after 
%                              each EM iteration.
%               LLcut       -- (1 x numIters) array; data log likelihood  
%                              after each EM iteration (where training data 
%                              was potentially 'cut' into trials of equal 
%                              length)
%                              Note: Entries will be NaN, on iterations
%                                    was LL was not computed, to save time.
%
%     Optional:
%     
%     freqLL    -- int; if known, specify how often LL was computed during
%                  model fitting (default: 1)
%     freqParam -- int; if known, specify how often GP parameters were
%                  stored during model fitting (default: 1)
%     units    -- string; units of time of binWidth (for labels)
%                 (default: '')
%
% Outputs:
%     None. (But creates a figure)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     12 Apr 2020 -- Initial full revision.
%     17 Apr 2020 -- Added 0-within-group dimension functionality
%     28 Jun 2020 -- Added 0-across-group dimension functionality
%     19 Feb 2023 -- Added spectral Gaussian compatibility

freqLL = 1;
freqParam = 1;
units = '';
assignopts(who,varargin);

% Format units for axis labels
if ~isempty(units)
    timeUnits = sprintf(' (%s)', units); 
    freqUnits = sprintf(' (1/%s)', units);
end

xDim_across = res.xDim_across;
xDim_within = res.xDim_within;
binWidth = res.binWidth;
numGroups = length(xDim_within);
numPlot = 2 + 2*(xDim_across > 0) + sum(xDim_within > 0);
if isequal(res.estParams.covType, 'sg')
    % Account for extra parameter in spectral Gaussian kernel
    numPlot = numPlot + (xDim_across > 0) + sum(xDim_within > 0); 
end

colors = generateColors(); % Get custom colors for plotting

figure;

% Log-likelihood curve
plotIdx = 1;
subplot(1,numPlot,plotIdx);
hold on;
LLcut = res.LLcut(~isnan(res.LLcut));
if freqLL > 2
    plot([1 2 freqLL.*(1:length(LLcut)-2)], LLcut, 'color', colors.grays{1}, 'linewidth', 1.5);
elseif freqLL == 2
    plot([1 freqLL.*(1:length(LLcut)-1)], LLcut, 'color', colors.grays{1}, 'linewidth', 1.5);
else
    plot(LLcut, 'color', colors.grays{1}, 'linewidth', 1.5);
end
xlabel('# Iterations');
ylabel('LL');

% Cumulative fitting time
plotIdx = plotIdx + 1;
subplot(1,numPlot,plotIdx);
hold on;
plot(cumsum(res.iterTime), 'color', colors.grays{1}, 'linewidth', 1.5);
xlabel('# Iterations');
ylabel('Cumulative time elapsed (sec)');

% Delay progress
if xDim_across > 0
    plotIdx = plotIdx + 1;
    subplot(1,numPlot,plotIdx);
    hold on;
    Delays = cat(3, res.D{:});
    for i = 1:xDim_across
        plot(freqParam.*(0:length(res.D)-1), binWidth.*squeeze(Delays(2,i,:)),...
             'linewidth', 1.5);
    end
    xlabel('# Iterations');
    ylabel(sprintf('Delay%s', timeUnits));

    % Across-group GP timescale progress
    Gams = cat(3, res.gams_across{:});
    plotIdx = plotIdx + 1;
    subplot(1,numPlot,plotIdx);
    hold on;
    for i = 1:xDim_across
        plot(freqParam.*(0:length(res.gams_across)-1), ...
             binWidth./sqrt(squeeze(Gams(1,i,:))), ...
             'linewidth', 1.5);
    end
    xlabel('# Iterations');
    ylabel(sprintf('Across-group GP timescales%s', timeUnits));
    
    if isequal(res.estParams.covType, 'sg')
        % Across-group GP center frequency progress
        Nus = cat(3, res.nus_across{:});
        plotIdx = plotIdx + 1;
        subplot(1,numPlot,plotIdx);
        hold on;
        for i = 1:xDim_across
            plot(freqParam.*(0:length(res.nus_across)-1), ...
                 squeeze(Nus(1,i,:))./binWidth, ...
                 'linewidth', 1.5);
        end
        xlabel('# Iterations');
        ylabel(sprintf('Across-group GP frequencies%s', freqUnits));
    end
end
    
% Within-group GP timescale progress
for groupIdx = 1:numGroups
    % Don't try to plot anything for groups with 0 within-group latents
    if xDim_within(groupIdx) > 0
        plotIdx = plotIdx + 1;
        subplot(1,numPlot,plotIdx);
        Gams = cat(3, res.gams_within{groupIdx}{:});
        hold on;
        for i = 1:xDim_within(groupIdx)
            plot(freqParam.*(0:length(res.gams_within{groupIdx})-1), ...
                 binWidth./sqrt(squeeze(Gams(1,i,:))), ...
                 'linewidth', 1.5);
        end
        xlabel('# Iterations');
        ylabel(sprintf('Within-group GP timescales, group %d%s', groupIdx, timeUnits))
    end
end

if isequal(res.estParams.covType, 'sg')
    % Within-group GP center frequency progress
    for groupIdx = 1:numGroups
        % Don't try to plot anything for groups with 0 within-group latents
        if xDim_within(groupIdx) > 0
            plotIdx = plotIdx + 1;
            subplot(1,numPlot,plotIdx);
            Nus = cat(3, res.nus_within{groupIdx}{:});
            hold on;
            for i = 1:xDim_within(groupIdx)
                plot(freqParam.*(0:length(res.nus_within{groupIdx})-1), ...
                     squeeze(Nus(1,i,:))./binWidth, ...
                     'linewidth', 1.5);
            end
            xlabel('# Iterations');
            ylabel(sprintf('Within-group GP frequencies, group %d%s', groupIdx, freqUnits))
        end
    end
end