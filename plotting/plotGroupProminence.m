function plotGroupProminence(prom, varargin)
%
% plotGroupProminence(prom,...)
%
% Description: Plot the bootstrapped "prominence" of each group of latents
%              among each group of observations in prom.
%
% Arguments:
%
%     Required:
%
%     prom -- (1 x numGroups) cell array; prom{i} contains a structure with
%             the bootstrapped prominence of each group of latent variables
%             for observation group i. See bootstrapProminence for details
%             on the format of these structures.
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
%     12 Jan 2021 -- Fixed issue with indexing that only arose for
%                    numGroups > 2.

metric    = 'LL';
extraOpts = assignopts(who, varargin);

colors = generateColors(); % Generate custom plotting colors

numGroups = length(prom);

% Determine number of subplots
figure;
for groupIdx = 1:numGroups
    currProm = prom{groupIdx};
    numAcross = length(currProm.across.dimGroups);
    numWithin = length(currProm.within.dimGroups{1});

    % Determine vertical axis limits for all subplots for this group
    upper = [currProm.across.(metric).upper currProm.within.(metric).upper{1}];
    lower = [currProm.across.(metric).lower currProm.within.(metric).lower{1}];
    upper = max(upper);
    lower = min([lower 0]); % Include 0, if all lower bounds are above 0
    range = upper - lower;
    vlim = ([lower-0.1*range upper+0.1*range]);

    % Across-group prominences
    subplot(numGroups,2,(groupIdx-1)*2+1);
    hold on;

    % Format subplot
    xlim([0 numAcross+1]);
    set(gca, 'XTick', 1:numAcross);
    xtcklbl = cell(1,numAcross);
    for i = 1:numAcross
        xtcklbl{i} = sprintf('[%s]', num2str(currProm.across.dimGroups{i})); 
    end
    set(gca, 'XTickLabel', xtcklbl);
    xlabel(sprintf('Across-group latents, group %d', groupIdx));
    ylim(vlim);
    ylabel(['Prominence (\Delta' sprintf('%s)', metric)]);

    % Plot data
    line([0 numAcross+1], [0 0], ...
         'Color', colors.grays{6}, ...
         'linestyle', '--', ...
         'linewidth', 1.5);
    errorbar(1:numAcross, ...
             currProm.across.(metric).raw, ...
             currProm.across.(metric).raw - currProm.across.(metric).lower, ...
             currProm.across.(metric).upper - currProm.across.(metric).raw, ...
             'o', ...
             'Color', colors.grays{1}, ...
             'linestyle', 'none', ...
             'linewidth', 1.5, ...
             'markerfacecolor', colors.grays{1});
    hold off;

    % Within-group prominences
    subplot(numGroups,2,(groupIdx-1)*2+2);
    hold on;

    % Format subplot
    xlim([0,numWithin+1]); 
    set(gca, 'XTick', 1:numWithin);
    xtcklbl = cell(1,numWithin);
    for i = 1:numWithin
        xtcklbl{i} = sprintf('[%s]', num2str(currProm.within.dimGroups{1}{i})); 
    end
    set(gca, 'XTickLabel', xtcklbl);
    xlabel(sprintf('Within-group latents, group %d', groupIdx)); 
    ylim(vlim);
    ylabel(['Prominence (\Delta' sprintf('%s)', metric)]);

    % Plot data
    line([0 numWithin+1], [0 0], ...
         'Color', colors.grays{6}, ...
         'linestyle', '--', ...
         'linewidth', 1.5);

    errorbar(1:numWithin, ...
             currProm.within.(metric).raw{1}, ...
             currProm.within.(metric).raw{1} - currProm.within.(metric).lower{1}, ...
             currProm.within.(metric).upper{1} - currProm.within.(metric).raw{1}, ...
             'o', ...
             'Color', colors.grays{1}, ...
             'linestyle', 'none', ...
             'linewidth', 1.5, ...
             'markerfacecolor', colors.grays{1});
    hold off;
end