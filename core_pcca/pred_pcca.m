function [R2, MSE] = pred_pcca(Ys, params)
%
% [R2, MSE] = pred_pcca(Ys, params)
%
% Description: Performs leave-group-out prediction using an existing pCCA
%              model.
%
% Arguments:
%     Ys     -- (1 x numGroups) cell array; list of data matrices 
%               {(y1Dim x N), (y2Dim x N), ...}
%     params -- learned pCCA parameters (structure with fields Cs, Rs, ds)
%
% Outputs:
%
%     R2.agg    -- float; aggregate leave-group-out R^2
%     R2.indiv  -- (1 x numGroups) array; R2.indiv(i) gives the R^2 value
%                  when predicting group i given the remaining groups
%     MSE.agg   -- float; aggregate leave-group-out mean-squared error
%     MSE.indiv -- (1 x numGroups) array; MSE.indiv(i) gives the MSE
%                  when predicting group i given the remaining groups
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     06 Aug 2023 -- Initial full revision.

numGroups = length(Ys);
N = size(Ys{1},2);
yDims = nan(1,numGroups);
Ys_pred = cell(1,numGroups);
for groupIdx = 1:numGroups
    yDims(groupIdx) = size(Ys{groupIdx},1);
    Ys_pred{groupIdx} = nan(size(Ys{groupIdx}));
end

% Perform leave-group-out prediction
R2.indiv = nan(1,numGroups);
MSE.indiv = nan(1,numGroups);
for groupIdx = 1:numGroups

    targetGroup = groupIdx; % Group to be left out
    sourceGroups = setdiff(1:numGroups, targetGroup); % Observed groups
    
    Ysource = vertcat(Ys{sourceGroups});
    Csource = vertcat(params.Cs{sourceGroups});
    Rsource = blkdiag(params.Rs{sourceGroups});
    dsource = vertcat(params.ds{sourceGroups});
    
    Ctarget = params.Cs{targetGroup};
    dtarget = params.ds{targetGroup};
    
    % Predict
    Ys_pred{groupIdx} = Ctarget * Csource' / (Csource * Csource' + Rsource) * (Ysource - repmat(dsource,1,N)) ...
                        + repmat(dtarget,1,N);
    
    % Compute performance metrics for this target group
    % MSE
    MSE.indiv(targetGroup) = immse(Ys_pred{targetGroup}, Ys{targetGroup});
    % R2
    RSS = sum( sum( ( Ys{targetGroup} - Ys_pred{targetGroup} ).^2, 1 ) );
    TSS = sum( sum( ( Ys{targetGroup} - repmat( mean(Ys{targetGroup},2), [1 size(Ys{targetGroup},2)] ) ).^2, 1 ) );
    R2.indiv(targetGroup) = 1 - RSS / TSS;

end

% Compute aggregate performance metrics
Ytrue = vertcat(Ys{:});
Ypred = vertcat(Ys_pred{:});
% MSE
MSE.agg = immse(Ypred, Ytrue);
% R2
RSS = sum( sum( ( Ytrue - Ypred ).^2, 1 ) );
TSS = sum( sum( ( Ytrue - repmat( mean(Ytrue,2), [1 size(Ytrue,2)] ) ).^2, 1 ) );
R2.agg = 1 - RSS / TSS;
