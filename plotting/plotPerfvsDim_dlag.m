function plotPerfvsDim_dlag(res,bestModel,plotAcross,varargin)
% 
% plotPerfvsDim_dlag(res,bestModel,plotAcross...)
%
% Description: Plot cross-validated performance metrics vs
%              latent dimensionality.
%
% Arguments:
%
%     Required:
%
%     res -- structure whose i-th entry has fields
%              xDim_across -- int; across-group latent dimensionality
%              xDim_within -- (1 x numGroups) array; within-group latent
%                             dimensionalities for each group
%              sumLL  -- float; cross-validated log-likelihood
%              LL     -- float; average cross-validated log-likelihood
%              LL_sem -- float; standard error of LL across CV folds
%              R2     -- (1 x 2) array; average cross-validated R^2 in each
%                        pairwise direction, for the pair in rGroups.
%              R2_sem -- (1 x 2) array; standard error of R2 across CV folds
%              R2orth -- (xDim_across x 2) array; average cross-validated
%                        R^2 error in each direction for reduced DLAG 
%                        predictions, for the pair in rGroups.
%              R2orth_sem -- (xDim x 2) array; standard error of R2orth
%                            across CV folds
%              MSE    -- (1 x 2) array; average cross-validated 
%                        mean-squared error in each pairwise direction, for 
%                        the pair in rGroups.
%              MSE_sem -- (1 x 2) array; standard error of MSE across CV 
%                         folds
%              MSEorth -- (xDim_across x 2) array; average cross-validated 
%                         mean-squared error in each direction for reduced 
%                         DLAG predictions, for the pair in rGroups.
%              MSEorth_sem -- (xDim_across x 2) array; standard error of
%                         MSE orth across CV folds
%              estParams -- model parameters estimated using all data
%              rGroups -- (1 x 2) array; the indexes of two groups 
%                         used to measure generalization performance
%                         via regression
%
%     bestModel -- int; index corresponding to the best model 
%                  (based on sumLL) in res (default:[])
%
%     plotAcross -- logical; set true to plot performance vs across-group
%                   dimensionality, while holding within-group
%                   dimensionalities fixed.
% 
%     Optional:
%
%     groupIdx  -- int; if plotAcross is false, then plot performance vs
%                  within-group dimensionality. groupIdx specifies
%                  which group to use. (default: 1)
%     fixWithin -- (1 x numGroups) array; Within-group dimensionalities at
%                  which to fix each group. If plotAcross is false, then 
%                  fixWithin(groupIdx) will be ignored, and that group's
%                  dimensionality will vary. 
%                  (default: res(bestModel).xDim_within)
%     fixAcross -- int; Fixed across-group dimensionality 
%                  (default: res(bestModel).xDim_across)
%     plotLL    -- logical; set true to plot data log likelihood 
%                  (default: true)
%     plotR2    -- logical; set true to plot R^2 (default: false)
%     plotMSE   -- logical; set true to plot MSE (default: false)
%     plotOrth  -- logical; set true to plot R^2 and MSE for the best
%                  performing DLAG model, orthonormalized. (default: false)
% 
% Outputs:
%     None. (But creates figures)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     09 Apr 2020 -- Initial full revision.

groupIdx = 1;
fixWithin = res(bestModel).xDim_within;
fixAcross = res(bestModel).xDim_across;
plotLL = true;
plotR2 = false;
plotMSE = false;
plotOrth = false;
assignopts(who, varargin);

numGroups = length(fixWithin);
numModels = length(res);
numRGroups = length(res(1).rGroups);   % Number of groups in a regression pair
colors = generateColors(); % Generate custom plotting colors

keepIdxs = [];
if plotAcross
    % Only plot models with within-group dimensionalities that match
    % fixWithin
    for modelIdx = 1:numModels
        if isequal(res(modelIdx).xDim_within, fixWithin)
            keepIdxs = [keepIdxs modelIdx]; 
        end
    end
    % Update bestModel to correspond to the kept models (if it remains at
    % all)
    bestModel = find(keepIdxs == bestModel);
    res = res(keepIdxs);
    numModels = length(res); % Update the number of remaining models
    xDim = [res.xDim_across];
    xlbl = {'Across-group dimensionality'; ...
            sprintf('(Within-group fixed at [%s])', num2str(fixWithin))};
else
    % Get the indices of all groups other than groupIdx
    otherIdxs = setdiff(1:numGroups, groupIdx);
    % Only plot models with across- and within-group dimensionalities that
    % match fixAcross and fixWithin
    xDim = [];
    for modelIdx = 1:numModels
        if isequal(res(modelIdx).xDim_within(otherIdxs), fixWithin(otherIdxs)) ...
                && isequal(res(modelIdx).xDim_across, fixAcross)
            
            keepIdxs = [keepIdxs modelIdx]; 
            xDim = [xDim res(modelIdx).xDim_within(groupIdx)];
        end
    end
    % Update bestModel to correspond to the kept models (if it remains at
    % all)
    bestModel = find(keepIdxs == bestModel);
    res = res(keepIdxs);
    numModels = length(res); % Update the number of remaining models
    xlbl = {sprintf('Within-group %d dimensionality', groupIdx); ...
            sprintf('(Across-group fixed at %d, Within-group fixed at [%s])',...
                   fixAcross, num2str(fixWithin(otherIdxs)))};
end



% Determine number of subplots, based on which metrics are to be plotted
numPlots = sum([plotLL plotR2 plotMSE]);

figure;
plotIdx = 0;

% Plot LL versus latent dimensionality
if plotLL
    plotIdx = plotIdx + 1;
    subplot(1,numPlots,plotIdx);
    hold on;
    xlabel(xlbl);
    ylabel('Cross-validated LL');
    
    sumLL = [res.sumLL];
    plot(xDim, sumLL, 'o-', 'Color', colors.grays{1}, ...
         'MarkerFaceColor', colors.grays{1}, 'linewidth', 1.5);
    if ~isempty(bestModel)
        legendEntries = plot(xDim(bestModel), sumLL(bestModel), 'p', ...
                             'color', colors.reds{4}, ...
                             'markerfacecolor', colors.reds{4}, ...
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
    xlabel(xlbl);
    ylabel('Cross-validated R^2');
    rcolors = {colors.grays{1}, colors.grays{6}};
    
    % Plot R^2 for prediction in both directions
    legendEntries = [];
    legendLabels = {};
    if ~plotOrth
        for rIdx = 1:numRGroups
            R2 = nan(1,numModels);
            R2_sem = nan(1,numModels);
            for modelIdx = 1:numModels
                R2(modelIdx) = res(modelIdx).R2(rIdx);
                R2_sem(modelIdx) = res(modelIdx).R2_sem(rIdx);
            end
            legendEntries(end+1) = errorbar(xDim, R2, R2_sem, 'o-', ...
                'linewidth', 1.5, 'Color', rcolors{rIdx}, ...
                'MarkerFaceColor', rcolors{rIdx});
            predGroup = setdiff(res(1).rGroups, res(1).rGroups(rIdx));
            legendLabels{end+1} = sprintf('Predicting group %01d', predGroup);
        end
    end
    
    if plotOrth && ~isempty(bestModel)
        xDim = 1:res(bestModel).xDim_across; % xDim potentially changes here
        for rIdx = 1:numRGroups
            R2orth = res(bestModel).R2orth(:,rIdx);
            R2orth_sem = res(bestModel).R2orth_sem(:,rIdx);
            legendEntries(end+1) = errorbar(xDim, R2orth, R2orth_sem, 'o-', ...
                'linewidth', 1.5, 'Color', rcolors{rIdx}, 'MarkerFaceColor', rcolors{rIdx});
            predGroup = setdiff(res(1).rGroups, res(1).rGroups(rIdx));
            legendLabels{end+1} = sprintf('Predicting group %01d', predGroup);
        end
    end
    legend(legendEntries, legendLabels, 'Location', 'southeast');
    hold off;
end

% Plot MSE versus latent dimensionality
if plotMSE
    plotIdx = plotIdx + 1;
    subplot(1,numPlots,plotIdx);
    hold on;
    xlabel(xlbl);
    ylabel('Cross-validated MSE');
    rcolors = {colors.grays{1}, colors.grays{6}};
    
    % Plot MSE for prediction in both directions
    legendEntries = [];
    legendLabels = {};
    if ~plotOrth
        for rIdx = 1:numRGroups
            MSE = nan(1,numModels);
            MSE_sem = nan(1,numModels);
            for modelIdx = 1:numModels
                MSE(modelIdx) = res(modelIdx).MSE(rIdx);
                MSE_sem(modelIdx) = res(modelIdx).MSE_sem(rIdx);
            end
            legendEntries(end+1) = errorbar(xDim, MSE, MSE_sem, 'o-', ...
                'linewidth', 1.5, 'Color', rcolors{rIdx}, ...
                'MarkerFaceColor', rcolors{rIdx});
            predGroup = setdiff(res(1).rGroups, res(1).rGroups(rIdx));
            legendLabels{end+1} = sprintf('Predicting group %01d', predGroup);
        end
    end
    
    if plotOrth && ~isempty(bestModel)
        xDim = 1:res(bestModel).xDim_across; % xDim potentially changes here
        for rIdx = 1:numRGroups
            MSEorth = res(bestModel).MSEorth(:,rIdx);
            MSEorth_sem = res(bestModel).MSEorth_sem(:,rIdx);
            legendEntries(end+1) = errorbar(xDim, MSEorth, MSEorth_sem, ...
                'o-', 'linewidth', 1.5, 'Color', rcolors{rIdx}, ...
                'MarkerFaceColor', rcolors{rIdx});
            predGroup = setdiff(res(1).rGroups, res(1).rGroups(rIdx));
            legendLabels{end+1} = sprintf('Predicting group %01d', predGroup);
        end
    end
    legend(legendEntries, legendLabels, 'Location', 'northeast');
    hold off;
end