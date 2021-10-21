function plotPerfvsDim_dlag(res,varargin)
% 
% plotPerfvsDim_dlag(res,...)
%
% Description: Plot cross-validated performance metrics vs
%              latent dimensionality.
%
% Arguments:
%
%     Required:
%
%     res -- structure whose i-th entry has (relevant) fields
%              xDim_across -- int; across-group latent dimensionality
%              xDim_within -- (1 x numGroups) array; within-group latent
%                             dimensionalities for each group
%              sumLL.joint -- float; cross-validated log-likelihood, 
%                             evaluated on all groups jointly
%              sumLL.indiv -- (1 x numGroups) array; cross-validated
%                             log-likelihood, evaluated on each group 
%                             individually
%              LL.joint    -- float; average cross-validated
%                             log-likelihood, evaluated on all groups
%                             jointly
%              LL.indiv    -- (1 x numGroups) array; average cross-validated
%                             log-likelihood, evaluated on each group 
%                             individually
%              LL_sem.joint -- float; standard error of LL.joint across CV
%                              folds
%              LL_sem.indiv -- (1 x numGroups) array; standard error of
%                              LL.indiv across CV folds
%              R2_reg.joint -- float; average cross-validated R^2,
%                              evaluated jointly across the pair in rGroups
%              R2_reg.indiv -- (1 x 2) array; average cross-validated R^2 
%                              in each pairwise direction, for the pair in 
%                              rGroups
%              R2_reg_sem.joint -- float; standard error of R2_reg.joint 
%                                  across CV folds
%              R2_reg_sem.indiv -- (1 x 2) array; standard error of 
%                                  R2_reg.indiv across CV folds
%              MSE_reg.indiv    -- (1 x 2) array; average cross-validated 
%                                  mean-squared error in each pairwise 
%                                  direction, for the pair in rGroups.
%              MSE_reg_sem.indiv -- (1 x 2) array; standard error of 
%                                   MSE_reg.indiv across CV folds
%              R2_denoise.joint -- float; average cross-validated R^2,
%                                  evaluated jointly from denoised 
%                                  reconstruction across all groups
%              R2_denoise.indiv -- (1 x numGroups) array; average 
%                                  cross-validated R^2, evaluated
%                                  for each group from denoised
%                                  reconstruction
%              R2_denoise_sem.joint -- float; standard error of
%                                      R2_denoise.joint across CV folds
%              R2_denoise_sem.indiv -- (1 x numGroups) array; standard error 
%                                      of R2_denoise.indiv across CV folds
%              MSE_denoise.joint -- float; average cross-validated 
%                                   mean-squared error, evaluated jointly
%                                   from denoised reconstruction across 
%                                   all groups
%              MSE_denoise.indiv -- (1 x numGroups) array; average 
%                                   cross-validated mean-squared error,
%                                   evaluated for each group from denoised
%                                   reconstruction
%              MSE_denoise_sem.joint -- float; standard error of
%                                      MSE_denoise.joint across CV folds
%              MSE_denoise_sem.indiv -- (1 x numGroups) array; standard error 
%                                      of MSE_denoise.indiv across CV folds
%              estParams -- model parameters estimated using all data
%              rGroups -- (1 x 2) array; the indexes of two groups 
%                         used to measure generalization performance
%                         via regression
% 
%     Optional:
%
%     verbose    -- logical; set true to plot additional performance metrics
%                   beyond cross-validated log-likelihood. These additional
%                   metrics may provide additional insight to more
%                   knowledgeable users, but may simply be confusing for
%                   the average user. (default: false)
%     plotAcross -- logical; set true to plot performance vs across-group
%                   dimensionality, while holding within-group
%                   dimensionalities fixed.
%     groupIdx  -- int; if plotAcross is false, then plot performance vs
%                  within-group dimensionality. groupIdx specifies
%                  which group to use. (default: 1)
%     fixWithin -- (1 x numGroups) array; Within-group dimensionalities at
%                  which to fix each group. If plotAcross is false, then 
%                  fixWithin(groupIdx) will be ignored, and that group's
%                  dimensionality will vary. 
%                  (default: res(1).xDim_within)
%     fixAcross -- int; Fixed across-group dimensionality 
%                  (default: res(1).xDim_across)
%     xDims_grid -- (numModels x numGroups+1) array; Each row specifies
%                   a particular model to be plotted, where 
%                   xDims_grid(i,:) = [xDim_across_i xDim_within_i]
%                   gives the across- and within-group dimensionalities
%                   for model i to be considered. (default: [])
%                   NOTE: If xDims_grid is specified, then it overrides all
%                         other optional arguments.
% 
% Outputs:
%     None. (But creates figures)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     09 Apr 2020 -- Initial full revision.
%     13 May 2020 -- 'bestModel' now always derived from models to be
%                    plotted. As a result, setting 'plotOrth' to true will 
%                    never result in blank R^2 and MSE plots, as before.
%     12 Jun 2020 -- Updated for expanded cross-validation metrics.
%                    Plotting of orthonormalized DLAG performance moved
%                    elsewhere.
%     16 Aug 2020 -- Added xDims_grid option.
%     20 Oct 2021 -- Minor plotting convention updates. Added verbose
%                    option.

verbose = false;
plotAcross = true;
groupIdx = 1;
fixWithin = res(1).xDim_within;
fixAcross = res(1).xDim_across;
xDims_grid = [];
assignopts(who, varargin);

numGroups = length(fixWithin);
numModels = length(res);
rGroups = res(1).rGroups;
colors = generateColors(); % Generate custom plotting colors

% ====================================================
% Determine which models' performance will be plotted
% ====================================================

keepIdxs = [];
if ~isempty(xDims_grid)
    % Only plot models with across- and within-group dimensionalities that 
    % match the rows of xDims_grid.
    xtklbl = {};
    for gridIdx = 1:size(xDims_grid,1)
        xDim_across = xDims_grid(gridIdx,1);
        xDim_within = xDims_grid(gridIdx,2:end);
        for modelIdx = 1:numModels
            if res(modelIdx).xDim_across == xDim_across && isequal(res(modelIdx).xDim_within, xDim_within)
                keepIdxs = [keepIdxs modelIdx];
                xtklbl{end+1} = sprintf('[%s]', num2str(xDims_grid(gridIdx,:)));
            end
        end
    end
    xlbl = {'Model candidate'};
    
    
else
    
    if plotAcross
        % Only plot models with within-group dimensionalities that match
        % fixWithin
        xtklbl = {};
        for modelIdx = 1:numModels
            if isequal(res(modelIdx).xDim_within, fixWithin)
                keepIdxs = [keepIdxs modelIdx]; 
                xtklbl{end+1} = num2str(res(modelIdx).xDim_across);
            end
        end
        xlbl = {'Across-group dimensionality'; ...
                sprintf('(Within-group fixed at [%s])', num2str(fixWithin))};
    else
        % Get the indices of all groups other than groupIdx
        otherIdxs = setdiff(1:numGroups, groupIdx);
        % Only plot models with across- and within-group dimensionalities that
        % match fixAcross and fixWithin
        xtklbl = {};
        for modelIdx = 1:numModels
            if isequal(res(modelIdx).xDim_within(otherIdxs), fixWithin(otherIdxs)) ...
                    && isequal(res(modelIdx).xDim_across, fixAcross)

                keepIdxs = [keepIdxs modelIdx]; 
                xtklbl{end+1} = num2str(res(modelIdx).xDim_within(groupIdx));
            end
        end
        xlbl = {sprintf('Within-group %d dimensionality', groupIdx); ...
                sprintf('(Across-group fixed at %d, Within-group fixed at [%s])',...
                       fixAcross, num2str(fixWithin(otherIdxs)))};
    end
    
end

% Throw out irrelevant models
res = res(keepIdxs);
numModels = length(res); % Update the number of remaining models

% Among the kept models, determine the best model
sumLL = nan(1,numModels);
for modelIdx = 1:numModels
    sumLL(modelIdx) = res(modelIdx).sumLL.joint;
end
[~, bestModel] = max(sumLL);

% Set up the figure and subplots
figure;
if verbose
    % Plot additional performance metrics
    numCol = 3;
    numRow = 1 + numGroups;
else
    % Plot only cross-validated log-likelihood, which is the recommended
    % metric for model selection
    numCol = 1;
    numRow = 1;
end

% ====================================
% Plot joint cross-validation metrics
% ====================================

% Plot LL versus latent dimensionality
plotIdx = 1;
subplot(numRow,numCol,plotIdx);
hold on;
xlabel(xlbl);
if verbose
    % Disambiguate from marginal likelihoods for each group
    ylabel('Cross-validated LL, joint');
else
    ylabel('Cross-validated LL');
end

sumLL = nan(1,numModels);
for modelIdx = 1:numModels
    sumLL(modelIdx) = res(modelIdx).sumLL.joint;
end
plot(1:numModels, sumLL, 'o-', 'Color', colors.grays{1}, ...
     'MarkerFaceColor', colors.grays{1}, 'linewidth', 1.5);

% Mark the best model among the plotted models.
legendEntries = plot(bestModel, sumLL(bestModel), 'p', ...
                     'color', colors.reds{4}, ...
                     'markerfacecolor', colors.reds{4}, ...
                     'markersize', 10);
legendLabels = 'best model';
legend(legendEntries, legendLabels, 'Location', 'southeast');
xticks(1:numModels);
xticklabels(xtklbl);
hold off;

% Additional performance metrics
if verbose
    % Plot R2 versus latent dimensionality
    plotIdx = plotIdx + 1;
    subplot(numRow,numCol,plotIdx);
    hold on;
    legendEntries = [];
    legendLabels = {'denoised', 'regression'};

    % Based on denoised prediction
    R2_denoise = nan(1,numModels);
    R2_denoise_sem = nan(1,numModels);
    for modelIdx = 1:numModels
        R2_denoise(modelIdx) = res(modelIdx).R2_denoise.joint;
        R2_denoise_sem(modelIdx) = res(modelIdx).R2_denoise_sem.joint;
    end
    legendEntries(end+1) = errorbar(1:numModels, R2_denoise, R2_denoise_sem, 'o-', ...
        'linewidth', 1.5, 'Color', colors.grays{1}, ...
        'MarkerFaceColor', colors.grays{1});

    % Based on pairwise regression
    R2_reg = nan(1,numModels);
    R2_reg_sem = nan(1,numModels);
    for modelIdx = 1:numModels
        R2_reg(modelIdx) = res(modelIdx).R2_reg.joint;
        R2_reg_sem(modelIdx) = res(modelIdx).R2_reg_sem.joint;
    end
    legendEntries(end+1) = errorbar(1:numModels, R2_reg, R2_reg_sem, 'o-', ...
        'linewidth', 1.5, 'Color', colors.reds{3}, ...
        'MarkerFaceColor', colors.reds{3});

    xlabel(xlbl);
    ylabel('Cross-validated R^2, joint');
    legend(legendEntries, legendLabels, 'Location', 'southeast');
    xticks(1:numModels);
    xticklabels(xtklbl);
    hold off;

    % Plot MSE versus latent dimensionality
    plotIdx = plotIdx + 1;
    subplot(numRow,numCol,plotIdx);
    hold on;
    legendEntries = [];
    legendLabels = {'denoised', 'regression'};

    % Based on denoised prediction
    MSE_denoise = nan(1,numModels);
    MSE_denoise_sem = nan(1,numModels);
    for modelIdx = 1:numModels
        MSE_denoise(modelIdx) = res(modelIdx).MSE_denoise.joint;
        MSE_denoise_sem(modelIdx) = res(modelIdx).MSE_denoise_sem.joint;
    end
    legendEntries(end+1) = errorbar(1:numModels, MSE_denoise, MSE_denoise_sem, 'o-', ...
        'linewidth', 1.5, 'Color', colors.grays{1}, ...
        'MarkerFaceColor', colors.grays{1});

    % Based on pairwise regression
    MSE_reg = nan(1,numModels);
    MSE_reg_sem = nan(1,numModels);
    for modelIdx = 1:numModels
        MSE_reg(modelIdx) = res(modelIdx).MSE_reg.joint;
        MSE_reg_sem(modelIdx) = res(modelIdx).MSE_reg_sem.joint;
    end
    legendEntries(end+1) = errorbar(1:numModels, MSE_reg, MSE_reg_sem, 'o-', ...
        'linewidth', 1.5, 'Color', colors.reds{3}, ...
        'MarkerFaceColor', colors.reds{3});

    xlabel(xlbl);
    ylabel('Cross-validated MSE, joint');
    legend(legendEntries, legendLabels, 'Location', 'northeast');
    xticks(1:numModels);
    xticklabels(xtklbl);
    hold off;

    % =========================================
    % Plot individual cross-validation metrics
    % =========================================
    for groupIdx = 1:numGroups

        % Plot LL versus latent dimensionality
        plotIdx = groupIdx*numCol + 1;
        subplot(numRow,numCol,plotIdx);
        hold on;
        xlabel(xlbl);
        ylabel(sprintf('Cross-validated LL, group %d', groupIdx));

        sumLL = nan(1,numModels);
        for modelIdx = 1:numModels
            sumLL(modelIdx) = res(modelIdx).sumLL.indiv(groupIdx);
        end
        plot(1:numModels, sumLL, 'o-', 'Color', colors.grays{1}, ...
             'MarkerFaceColor', colors.grays{1}, 'linewidth', 1.5);
        xticks(1:numModels);
        xticklabels(xtklbl);
        hold off;

        % Plot R2 versus latent dimensionality
        plotIdx = plotIdx + 1;
        subplot(numRow,numCol,plotIdx);
        hold on;
        legendEntries = [];
        legendLabels = {};

        % Based on denoised prediction
        R2_denoise = nan(1,numModels);
        R2_denoise_sem = nan(1,numModels);
        for modelIdx = 1:numModels
            R2_denoise(modelIdx) = res(modelIdx).R2_denoise.indiv(groupIdx); 
            R2_denoise_sem(modelIdx) = res(modelIdx).R2_denoise_sem.indiv(groupIdx);
        end
        legendEntries(end+1) = errorbar(1:numModels, R2_denoise, R2_denoise_sem, 'o-', ...
            'linewidth', 1.5, 'Color', colors.grays{1}, ...
            'MarkerFaceColor', colors.grays{1});
        legendLabels{end+1} = 'denoised';

        % Based on pairwise regression
        if ismember(groupIdx,rGroups)
            % Get predictions conditioned on the other group in the pair
            predGroup = setdiff(rGroups, groupIdx);
            predIdx = find(rGroups == predGroup);
            R2_reg = nan(1,numModels);
            R2_reg_sem = nan(1,numModels);
            for modelIdx = 1:numModels
                R2_reg(modelIdx) = res(modelIdx).R2_reg.indiv(predIdx); 
                R2_reg_sem(modelIdx) = res(modelIdx).R2_reg_sem.indiv(predIdx);
            end
            legendEntries(end+1) = errorbar(1:numModels, R2_reg, R2_reg_sem, 'o-', ...
                'linewidth', 1.5, 'Color', colors.reds{3}, ...
                'MarkerFaceColor', colors.reds{3});
            legendLabels{end+1} = 'regression';
        end
        xlabel(xlbl);
        ylabel(sprintf('Cross-validated R^2, group %d', groupIdx));
        legend(legendEntries, legendLabels, 'Location', 'southeast');
        xticks(1:numModels);
        xticklabels(xtklbl);
        hold off;

        % Plot MSE versus latent dimensionality
        plotIdx = plotIdx + 1;
        subplot(numRow,numCol,plotIdx);
        hold on;
        legendEntries = [];
        legendLabels = {};

        % Based on denoised prediction
        MSE_denoise = nan(1,numModels);
        MSE_denoise_sem = nan(1,numModels);
        for modelIdx = 1:numModels
            MSE_denoise(modelIdx) = res(modelIdx).MSE_denoise.indiv(groupIdx); 
            MSE_denoise_sem(modelIdx) = res(modelIdx).MSE_denoise_sem.indiv(groupIdx);
        end
        legendEntries(end+1) = errorbar(1:numModels, MSE_denoise, MSE_denoise_sem, 'o-', ...
            'linewidth', 1.5, 'Color', colors.grays{1}, ...
            'MarkerFaceColor', colors.grays{1});
        legendLabels{end+1} = 'denoised';

        % Based on pairwise regression
        if ismember(groupIdx,rGroups)
            % Get predictions conditioned on the other group in the pair
            predGroup = setdiff(rGroups, groupIdx);
            predIdx = find(rGroups == predGroup);
            MSE_reg = nan(1,numModels);
            MSE_reg_sem = nan(1,numModels);
            for modelIdx = 1:numModels
                MSE_reg(modelIdx) = res(modelIdx).MSE_reg.indiv(predIdx); 
                MSE_reg_sem(modelIdx) = res(modelIdx).MSE_reg_sem.indiv(predIdx);
            end
            legendEntries(end+1) = errorbar(1:numModels, MSE_reg, MSE_reg_sem, 'o-', ...
                'linewidth', 1.5, 'Color', colors.reds{3}, ...
                'MarkerFaceColor', colors.reds{3});
            legendLabels{end+1} = 'regression';
        end
        xlabel(xlbl);
        ylabel(sprintf('Cross-validated MSE, group %d', groupIdx));
        legend(legendEntries, legendLabels, 'Location', 'northeast');
        xticks(1:numModels);
        xticklabels(xtklbl);
        hold off;
    end
    
end