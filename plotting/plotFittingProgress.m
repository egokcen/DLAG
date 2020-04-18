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

freqLL = 1;
freqParam = 1;
units = '';
assignopts(who,varargin);

% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

xDim_across = res.xDim_across;
xDim_within = res.xDim_within;
binWidth = res.binWidth;
numGroups = length(xDim_within);
numPlot = 4 + sum(xDim_within > 0);

colors = generateColors(); % Get custom colors for plotting

figure;

% Log-likelihood curve
subplot(1,numPlot,1);
hold on;
LLcut = res.LLcut(~isnan(res.LLcut));
plot([1 2 freqLL.*(1:length(LLcut)-2)], LLcut, 'color', colors.grays{1}, 'linewidth', 1.5);
xlabel('# Iterations');
ylabel('LL');

% Cumulative fitting time
subplot(1,numPlot,2);
hold on;
plot(cumsum(res.iterTime), 'color', colors.grays{1}, 'linewidth', 1.5);
xlabel('# Iterations');
ylabel('Cumulative time elapsed (sec)');

% Delay progress
subplot(1,numPlot,3);
hold on;
Delays = cat(3, res.D{:});
for i = 1:xDim_across
    plot(freqParam.*(0:length(res.D)-1), binWidth.*squeeze(Delays(2,i,:)),...
         'linewidth', 1.5);
end
xlabel('# Iterations');
ylabel(sprintf('Delay%s', units));

% Across-group GP timescale progress
Gams = cat(3, res.gams_across{:});
subplot(1,numPlot,4);
hold on;
for i = 1:xDim_across
    plot(freqParam.*(0:length(res.gams_across)-1), ...
         binWidth./sqrt(squeeze(Gams(1,i,:))), ...
         'linewidth', 1.5);
end
xlabel('# Iterations');
ylabel(sprintf('Across-group GP timescales%s', units));

% Within-group GP timescale progress
plotIdx = 4;
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
        ylabel(sprintf('Within-group GP timescales, group %d%s', groupIdx, units))
    end
end