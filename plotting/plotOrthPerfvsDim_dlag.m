function plotOrthPerfvsDim_dlag(res)
% 
% plotOrthPerfvsDim_dlag(res)
%
% Description: Plot cross-validated performance metrics vs
%              latent dimensionality for a single orthonormalized DLAG
%              model.
%
% Arguments:
%
%     res -- structure with the following (relevant) fields
%              xDim_across -- int; across-group latent dimensionality
%              xDim_within -- (1 x numGroups) array; within-group latent
%                             dimensionalities for each group
%              R2orth_reg.indiv -- (xDim_across x 2) array; average 
%                                  cross-validated R^2 error in each 
%                                  direction for orthonormalized DLAG 
%                                  predictions, for the pair in rGroups.
%              R2orth_reg_sem.indiv -- (xDim_across x 2) array; standard error 
%                                      of R2orth_reg.indiv across CV folds
%              MSEorth_reg.indiv -- (xDim_across x 2) array; average 
%                                   cross-validated mean-squared error in 
%                                   each direction for orthonormalized DLAG 
%                                   predictions, for the pair in rGroups.
%              MSEorth_reg_sem.indiv -- (xDim_across x 2) array; standard 
%                                       error of MSEorth_reg.indiv across
%                                       CV folds
%              R2orth_denoise.indiv -- (1 x numGroups) cell array;
%                                      average cross-validated R^2 for the 
%                                      orthonormalized DLAG model, evaluated
%                                      for each group from denoised
%                                      reconstruction
%              R2orth_denoise_sem.indiv -- (1 x numGroups) cell array; 
%                                      standard error of 
%                                      R2orth_denoise.indiv across CV folds
%              MSEorth_denoise.indiv -- (1 x numGroups) cell array;
%                                      average cross-validated mean-squared 
%                                      error for the orthonormalized DLAG 
%                                      model, evaluated for each group from
%                                      denoised reconstruction
%              MSEorth_denoise_sem.indiv -- (1 x numGroups) cell array; 
%                                      standard error of 
%                                      MSEorth_denoise.indiv across CV folds
%              estParams -- model parameters estimated using all data
%              rGroups -- (1 x 2) array; the indexes of two groups 
%                         used to measure generalization performance
%                         via regression
% 
% Outputs:
%     None. (But creates figures)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     13 Jun 2020 -- Initial full revision.
%     28 Jun 2020 -- Added 0-across-group dimension functionality

xDim_across = res.xDim_across;
xDim_within = res.xDim_within;
numGroups = length(xDim_within);
rGroups = res.rGroups;
colors = generateColors(); % Generate custom plotting colors

% Set up the figure and subplots
figure;
numCol = 2;
numRow = numGroups;
xlbl = 'Latent dimensionality';

% =========================================
% Plot individual cross-validation metrics
% =========================================
for groupIdx = 1:numGroups
    
    % Plot R2 versus latent dimensionality
    plotIdx = (groupIdx-1)*numCol + 1;
    subplot(numRow,numCol,plotIdx);
    hold on;
    legendEntries = [];
    legendLabels = {};
    
    % Based on denoised prediction
    % Including all types of latents
    legendEntries(end+1) = errorbar(0:length(res.R2orth_denoise.indiv{groupIdx})-1, ...
        res.R2orth_denoise.indiv{groupIdx}, res.R2orth_denoise_sem.indiv{groupIdx}, 'o-', ...
        'linewidth', 1.5, 'Color', colors.grays{1}, ...
        'MarkerFaceColor', colors.grays{1});
    legendLabels{end+1} = 'denoised';

    % Based on pairwise regression
    if ismember(groupIdx,rGroups)
        % Get predictions conditioned on the other group in the pair
        predGroup = setdiff(rGroups, groupIdx);
        predIdx = find(rGroups == predGroup);
        legendEntries(end+1) = errorbar(0:length(res.R2orth_reg.indiv(:,predIdx))-1, ...
            res.R2orth_reg.indiv(:,predIdx), res.R2orth_reg_sem.indiv(:,predIdx), 'o-', ...
            'linewidth', 1.5, 'Color', colors.reds{3}, ...
            'MarkerFaceColor', colors.reds{3});
        legendLabels{end+1} = 'regression';
    end
    xlabel(xlbl);
    ylabel(sprintf('Cross-validated R^2, group %d', groupIdx));
    legend(legendEntries, legendLabels, 'Location', 'southeast');
    hold off;

    % Plot MSE versus latent dimensionality
    plotIdx = plotIdx + 1;
    subplot(numRow,numCol,plotIdx);
    hold on;
    legendEntries = [];
    legendLabels = {};

    % Based on denoised prediction
    % Including all types of latents
    legendEntries(end+1) = errorbar(0:length(res.MSEorth_denoise.indiv{groupIdx})-1, ...
        res.MSEorth_denoise.indiv{groupIdx}, res.MSEorth_denoise_sem.indiv{groupIdx}, 'o-', ...
        'linewidth', 1.5, 'Color', colors.grays{1}, ...
        'MarkerFaceColor', colors.grays{1});
    legendLabels{end+1} = 'denoised';

    % Based on pairwise regression
    if ismember(groupIdx,rGroups)
        % Get predictions conditioned on the other group in the pair
        predGroup = setdiff(rGroups, groupIdx);
        predIdx = find(rGroups == predGroup);
        legendEntries(end+1) = errorbar(0:length(res.MSEorth_reg.indiv(:,predIdx))-1, ...
            res.MSEorth_reg.indiv(:,predIdx), res.MSEorth_reg_sem.indiv(:,predIdx), 'o-', ...
            'linewidth', 1.5, 'Color', colors.reds{3}, ...
            'MarkerFaceColor', colors.reds{3});
        legendLabels{end+1} = 'regression';
    end
    xlabel(xlbl);
    ylabel(sprintf('Cross-validated MSE, group %d', groupIdx));
    legend(legendEntries, legendLabels, 'Location', 'northeast');
    hold off;
    
end