function plotPerfvsDim_pcca(res,varargin)
% 
% plotPerfvsDim_pcca(res,...)
%
% Description: Plot cross-validated performance metrics vs
%              latent dimensionality.
%
% Arguments:
%
%     Required:
%
%     res -- structure whose i-th entry has fields
%              xDim   -- int; latent dimensionality
%              sumLL  -- float; cross-validated log-likelihood
%              LL     -- float; average cross-validated log-likelihood
%              LL_sem -- float; standard error of LL across CV folds
%              R2     -- (1 x 2) array; average cross-validated R^2 in each
%                        pairwise direction, for the pair in rGroups.
%              R2_sem -- (1 x 2) array; standard error of R2 across CV folds
%              MSE    -- (1 x 2) array; average cross-validated 
%                        mean-squared error in each pairwise direction, for 
%                        the pair in rGroups.
%              MSE_sem -- (1 x 2) array; standard error of MSE across CV 
%                         folds
%              estParams -- model parameters estimated using all data
%              rGroups -- (1 x 2) array; the indexes of two groups 
%                         used to measure generalization performance
%                         via regression
% 
%     Optional:
%
%     bestModel -- int; index corresponding to the best model 
%                  (based on sumLL) in res (default:[])
%     plotLL    -- logical; set true to plot data log likelihood 
%                  (default: true)
%     plotR2    -- logical; set true to plot R^2 (default: false)
%     plotMSE   -- logical; set true to plot MSE (default: false)
% 
% Outputs:
%     None. (But creates figures)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     09 Apr 2020 -- Initial full revision.

bestModel = [];
plotLL = true;
plotR2 = false;
plotMSE = false;
assignopts(who, varargin);

xDim = [res.xDim];
numModels = length(res);
numRGroups = length(res(1).rGroups);   % Number of groups in a regression pair
colors = generateColors(); % Generate custom plotting colors

% Determine number of subplots, based on which metrics are to be plotted
numPlots = sum([plotLL plotR2 plotMSE]);

figure;
plotIdx = 0;

% Plot LL versus latent dimensionality
if plotLL
    plotIdx = plotIdx + 1;
    subplot(1,numPlots,plotIdx);
    hold on;
    xlabel('Latent dimensionality');
    ylabel('Cross-validated LL');
    
    sumLL = [res.sumLL];
    plot(xDim, sumLL, 'o-', 'color', colors.grays{1}, ...
         'MarkerFaceColor', colors.grays{1}, 'linewidth', 1.5);
    if ~isempty(bestModel)
        legendEntries = plot(xDim(bestModel), sumLL(bestModel), 'p', ...
             'color', colors.reds{4}, ...
             'markerfacecolor', colors.reds{4},...
             'markersize', 10);
         legendLabels = 'best model';
         legend(legendEntries, legendLabels, 'Location', 'southeast');
    end
    hold off;
end

% Plot R2 versus latent dimensionality
if plotR2
    plotIdx = plotIdx + 1;
    subplot(1,numPlots,plotIdx);
    hold on;
    xlabel('Latent dimensionality');
    ylabel('Cross-validated R^2');
    rcolors = {colors.grays{1}, colors.grays{6}};
    
    % Plot R^2 for prediction in both directions
    legendEntries = [];
    legendLabels = {};
    for rIdx = 1:numRGroups
        R2 = nan(1,numModels);
        R2_sem = nan(1,numModels);
        for modelIdx = 1:numModels
            R2(modelIdx) = res(modelIdx).R2(rIdx);
            R2_sem(modelIdx) = res(modelIdx).R2_sem(rIdx);
        end
        legendEntries(end+1) = errorbar(xDim, R2, R2_sem, 'o-', ...
            'color', rcolors{rIdx}, 'MarkerFaceColor', rcolors{rIdx}, ...
            'linewidth', 1.5);
        predGroup = setdiff(res(1).rGroups, res(1).rGroups(rIdx));
            legendLabels{end+1} = sprintf('Predicting group %01d', predGroup);
    end
    legend(legendEntries, legendLabels, 'Location', 'southeast');
    hold off;
    
end

% Plot MSE versus latent dimensionality
if plotMSE
    plotIdx = plotIdx + 1;
    subplot(1,numPlots,plotIdx);
    hold on;
    xlabel('Latent dimensionality');
    ylabel('Cross-validated MSE');
    rcolors = {colors.grays{1}, colors.grays{6}};
    
    % Plot MSE for prediction in both directions
    legendEntries = [];
    legendLabels = {};
    for rIdx = 1:numRGroups
        MSE = nan(1,numModels);
        MSE_sem = nan(1,numModels);
        for modelIdx = 1:numModels
            MSE(modelIdx) = res(modelIdx).MSE(rIdx);
            MSE_sem(modelIdx) = res(modelIdx).MSE_sem(rIdx);
        end
        legendEntries(end+1) = errorbar(xDim, MSE, MSE_sem, 'o-', ...
            'color', rcolors{rIdx}, 'MarkerFaceColor', rcolors{rIdx}, ...
            'linewidth', 1.5);
        predGroup = setdiff(res(1).rGroups, res(1).rGroups(rIdx));
            legendLabels{end+1} = sprintf('Predicting group %01d', predGroup);
    end
    legend(legendEntries, legendLabels, 'Location', 'northeast');
    hold off;
end