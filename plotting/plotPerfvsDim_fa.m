function plotPerfvsDim_fa(res,varargin)
% 
% plotPerfvsDim_fa(res,...)
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
%              estParams -- model parameters estimated using all data
% 
%     Optional:
%
%     bestModels -- (1 x numGroups) array; indices corresponding to the 
%                   best models (based on sumLL) fit to each group (default:[])
% 
% Outputs:
%     None. (But creates figures)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

bestModels = [];
assignopts(who, varargin);
numGroups = length(res);

figure;

for groupIdx = 1:numGroups
    xDim = [res{groupIdx}.xDim];
    numModels = length(res{groupIdx});
    colors = generateColors(); % Generate custom plotting colors

    % Plot LL versus latent dimensionality
    subplot(1,numGroups,groupIdx);
    hold on;
    xlabel(sprintf('Latent dimensionality, group %d', groupIdx));
    ylabel('Cross-validated LL');

    sumLL = [res{groupIdx}.sumLL];
    plot(xDim, sumLL, 'o-', 'color', colors.grays{1}, ...
         'MarkerFaceColor', colors.grays{1}, 'linewidth', 1.5);
    if ~isempty(bestModels)
        legendEntries = plot(xDim(bestModels(groupIdx)), sumLL(bestModels(groupIdx)), 'p', ...
             'color', colors.reds{4}, ...
             'markerfacecolor', colors.reds{4},...
             'markersize', 10);
         legendLabels = 'best model';
         legend(legendEntries, legendLabels, 'Location', 'southeast');
    end
    hold off;
end