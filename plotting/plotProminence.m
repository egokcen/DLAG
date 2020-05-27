function plotProminence(prom, varargin)
%
% plotProminence(prom,...)
%
% Description: Plot the bootstrapped "prominence" of each group of latents
%              in prom.
%
% Arguments:
%
%     Required:
%
%     prom -- structure containing the following fields:
%             across -- structure with across-group prominence. Has fields
%                 LL.raw    -- (1 x numDimGroups) array; prominence of
%                              dimGroups_across(i) evaluated on raw data, 
%                              measured by decrease in log-likelihood 
%                              relative to the full model.
%                 LL.upper  -- (1 x numDimGroups) array; upper bound of
%                              bootstrap confidence interval
%                 LL.lower  -- (1 x numDimGroups) array; lower bound of
%                              bootstrap confidence interval
%                 VE.raw    -- (1 x numDimGroups) array;prominence of
%                              dimGroups_across(i) evaluated on raw data, 
%                              measured by normalized decrease in variance 
%                              explained relative to the full model.
%                 VE.upper  -- (1 x numDimGroups) array; upper bound of
%                              bootstrap confidence interval
%                 VE.lower  -- (1 x numDimGroups) array; lower bound of
%                              bootstrap confidence interval
%                 dimGroups -- (1x numDimGroups) cell array; Exactly
%                              dimGroups_across. Entries correspond to
%                              above fields.
%             within -- structure with within-group prominence. Has fields
%                 LL.raw    -- (1 x numGroups) cell array; prominence of 
%                              dimGroups_within{i} evaluated on raw data, 
%                              measured by decrease in log-likelihood 
%                              relative to the full model.
%                 LL.upper  -- (1 x numGroups) cell array; upper bound of
%                              bootstrap confidence interval
%                 LL.lower  -- (1 x numGroups) cell array; lower bound of
%                              bootstrap confidence interval
%                 VE.raw    -- (1 x numGroups) cell array; prominence of
%                              dimGroups_within{i} evaluated on raw data, 
%                              measured by normalized decrease in variance
%                              explained relative to the full model.
%                 VE.upper  -- (1 x numGroups) cell array; upper bound of
%                              bootstrap confidence interval
%                 VE.lower  -- (1 x numGroups) cell array; lower bound of
%                              bootstrap confidence interval
%                 dimGroups -- (1x numGroups) cell array; Exactly
%                              dimGroups_within. Entries correspond to
%                              above fields.
%
%     Optional:
%
%     metric -- string; Specify which measure of prominence to plot: 'LL'
%               or 'VE' (default: 'LL')
%
% Outputs:
%
%     None. (But creates figures)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 May 2020 -- Initial full revision.

metric    = 'LL';
extraOpts = assignopts(who, varargin);

numAcross = length(prom.across.dimGroups);
numGroups = length(prom.within.(metric).raw);
numWithin = nan(1,numGroups);
for groupIdx = 1:numGroups
    numWithin(groupIdx) = length(prom.within.dimGroups{groupIdx});
end

colors = generateColors(); % Generate custom plotting colors

% Determine number of subplots
numPlot = numGroups + 1;

% Determine vertical axis limits for all subplots
upper = prom.across.(metric).upper;
lower = prom.across.(metric).lower;
for groupIdx = 1:numGroups
    upper = [upper prom.within.(metric).upper{groupIdx}];
    lower = [lower prom.within.(metric).lower{groupIdx}];
end
upper = max(upper);
lower = min([lower 0]); % Include 0, if all lower bounds are above 0
range = upper - lower;
vlim = ([lower-0.1*range upper+0.1*range]);

figure;
plotIdx = 0;

% Across-group prominences
plotIdx = plotIdx + 1;
subplot(1,numPlot,plotIdx);
hold on;

% Format subplot
xlim([0 numAcross+1]);
set(gca, 'XTick', 1:numAcross);
xtcklbl = cell(1,numAcross);
for i = 1:numAcross
    xtcklbl{i} = sprintf('[%s]', num2str(prom.across.dimGroups{i})); 
end
set(gca, 'XTickLabel', xtcklbl);
xlabel('Across-group latents');
ylim(vlim);
ylabel(['Prominence (\Delta' sprintf('%s)', metric)]);

% Plot data
line([0 numAcross+1], [0 0], ...
     'Color', colors.grays{6}, ...
     'linestyle', '--', ...
     'linewidth', 1.5);
errorbar(1:numAcross, ...
         prom.across.(metric).raw, ...
         prom.across.(metric).raw - prom.across.(metric).lower, ...
         prom.across.(metric).upper - prom.across.(metric).raw, ...
         'o', ...
         'Color', colors.grays{1}, ...
         'linestyle', 'none', ...
         'linewidth', 1.5, ...
         'markerfacecolor', colors.grays{1});
hold off;

% Within-group prominences
for groupIdx = 1:numGroups
    plotIdx = plotIdx + 1;
    subplot(1,numPlot,plotIdx);
    hold on;
    
    % Format subplot
    xlim([0,numWithin(groupIdx)+1]); 
    set(gca, 'XTick', 1:numWithin(groupIdx));
    xtcklbl = cell(1,numWithin(groupIdx));
    for i = 1:numWithin(groupIdx)
        xtcklbl{i} = sprintf('[%s]', num2str(prom.within.dimGroups{groupIdx}{i})); 
    end
    set(gca, 'XTickLabel', xtcklbl);
    xlabel(sprintf('Within-group latents, group %d', groupIdx)); 
    ylim(vlim);
    ylabel(['Prominence (\Delta' sprintf('%s)', metric)]);
    
    % Plot data
    line([0 numWithin(groupIdx)+1], [0 0], ...
         'Color', colors.grays{6}, ...
         'linestyle', '--', ...
         'linewidth', 1.5);
 
    errorbar(1:numWithin(groupIdx), ...
             prom.within.(metric).raw{groupIdx}, ...
             prom.within.(metric).raw{groupIdx} - prom.within.(metric).lower{groupIdx}, ...
             prom.within.(metric).upper{groupIdx} - prom.within.(metric).raw{groupIdx}, ...
             'o', ...
             'Color', colors.grays{1}, ...
             'linestyle', 'none', ...
             'linewidth', 1.5, ...
             'markerfacecolor', colors.grays{1});
    hold off;
end