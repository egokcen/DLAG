function bestModel = getNumAcrossDim_dlag(res,modelList,varargin)
% 
% dim = getNumAcrossDim_dlag(res,modelList...)
%
% Description: Determine the optimal number of across-group dimensions 
%              required to explain interactions across groups, among the
%              model candidates specified in model list.
%
% Arguments:
%
%     Required:
%
%     res -- structure whose i-th entry has the following (relevant) fields
%              xDim_across -- int; across-group latent dimensionality
%              xDim_within -- (1 x numGroups) array; within-group latent
%                             dimensionalities for each group
%              R2_reg.joint -- float; average cross-validated R^2,
%                              evaluated jointly across the pair in rGroups
%              R2_reg.indiv -- (1 x 2) array; average cross-validated R^2 
%                              in each pairwise direction, for the pair in 
%                              rGroups
%              R2_reg_sem.joint -- float; standard error of R2_reg.joint 
%                                  across CV folds
%              R2_reg_sem.indiv -- (1 x 2) array; standard error of 
%                                  R2_reg.indiv across CV folds
%              MSE_reg.joint -- float; average cross-validated mean-squared
%                               error, evaluated jointly across the pair in
%                               rGroups
%              MSE_reg.indiv     -- (1 x 2) array; average cross-validated 
%                                   mean-squared error in each pairwise 
%                                   direction, for the pair in rGroups.
%              MSE_reg_sem.joint -- float; standard error of MSE_reg.joint 
%                                   across CV folds
%              MSE_reg_sem.indiv -- (1 x 2) array; standard error of 
%                                   MSE_reg.indiv across CV folds
%              estParams -- model parameters estimated using all data
%              rGroups -- (1 x 2) array; the indexes of two groups 
%                         used to measure generalization performance
%                         via regression
%     modelList -- (numModels x numGroups+1) array; modelList(i,:) gives
%                  the desired number of across- and within-group
%                  dimensions for model i. The first column of
%                  modelList(i,:) specifies the number of across-group
%                  dimensions. The remaining columns specify the number
%                  of within-group dimensions for each group. For
%                  example,
%                      modelList = [1 2 3; 4 5 6]
%                  specifies model 1 with xDim_across = 1,
%                  xDim_within = [2 3]; model 2 has xDim_across = 4, 
%                  xDim_within = [5 6].
%
%     Optional:
%
%     metric -- string; choose which cross-validated performance metric to
%               use for selection ('MSE', 'R2') (default: 'R2')
%     joint  -- logical; if true, select the optimal model according to a 
%               joint metric. If false, choose the optimal model according
%               to the minimum value among all individual metrics 
%               (default: false)
% 
% Outputs:
%     bestModel -- structure corresponding to the optimal model in res.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Aug 2020 -- Initial full revision.

metric = 'R2';
joint = false;
assignopts(who, varargin);

numModels = length(modelList);
numRGroups = 2; % Individual regression metrics are pairwise: one prediction in each direction

% Collect all dimensionalities contained in res
xDims = [];
for modelIdx = 1:length(res)
    xDims(modelIdx,:) = [res(modelIdx).xDim_across res(modelIdx).xDim_within]; 
end

% Consider only the models specified in modelList
keptModels = find(ismember(xDims, modelList, 'rows'));
resKept = res(keptModels);

if isequal(metric,'MSE')
    
    if joint
        % Collect performance metrics across models
        MSE_joint = nan(numModels,1);
        MSE_joint_sem = nan(numModels,1);
        for modelIdx = 1:numModels
            MSE_joint(modelIdx) = resKept(modelIdx).MSE_reg.joint;
            MSE_joint_sem(modelIdx) = resKept(modelIdx).MSE_reg_sem.joint;
        end

        % Determine the number of across-group dimensions
        [bestMSE, bestIdx] = min(MSE_joint);
        best_sem = MSE_joint_sem(bestIdx);
        % NOTE: max used in this manner selects the first instance that
        %       satisfies the input condition.
        [~, bestIdx] = max(MSE_joint <= bestMSE + best_sem); 
        bestModel = resKept(bestIdx);
        
    else
        % Collect performance metrics across models
        MSE_indiv = nan(numModels,numRGroups);
        MSE_indiv_sem = nan(numModels,numRGroups);
        for modelIdx = 1:numModels
            MSE_indiv(modelIdx,:) = resKept(modelIdx).MSE_reg.indiv;
            MSE_indiv_sem(modelIdx,:) = resKept(modelIdx).MSE_reg_sem.indiv;
        end

        % Determine the number of across-group dimensions required to 
        % predict the other group.
        bestIdxs = nan(1,numRGroups); % Indexes of best models in resKept.
        bestAcross = nan(1,numRGroups); % Across-group dimensionalities of best models
        for groupIdx = 1:numRGroups
            [bestMSE, bestIdx] = min(MSE_indiv(:,groupIdx));
            best_sem = MSE_indiv_sem(bestIdx,groupIdx);
            % NOTE: max used in this manner selects the first instance that
            %       satisfies the input condition.
            [~, bestIdxs(groupIdx)] = max(MSE_indiv(:,groupIdx) <= bestMSE + best_sem); 
            bestAcross(groupIdx) = resKept(bestIdxs(groupIdx)).xDim_across;
        end

        % Take the model with the minimum of the optimal values in each direction
        [~, bestGroupIdx] = min(bestAcross);   % Gives index into bestIdxs.
        bestModelIdx = bestIdxs(bestGroupIdx); % Gives index into resKept.
        bestModel = resKept(bestModelIdx);
    end

elseif isequal(metric,'R2')
    if joint
        % Collect performance metrics across models
        R2_joint = nan(numModels,1);
        R2_joint_sem = nan(numModels,1);
        for modelIdx = 1:numModels
            R2_joint(modelIdx) = resKept(modelIdx).R2_reg.joint;
            R2_joint_sem(modelIdx) = resKept(modelIdx).R2_reg_sem.joint;
        end

        % Determine the number of across-group dimensions
        [bestR2, bestIdx] = max(R2_joint);
        best_sem = R2_joint_sem(bestIdx);
        % NOTE: max used in this manner selects the first instance that
        %       satisfies the input condition.
        [~, bestIdx] = max(R2_joint >= bestR2 - best_sem); 
        bestModel = resKept(bestIdx);
        
    else
        % Collect performance metrics across models
        R2_indiv = nan(numModels,numRGroups);
        R2_indiv_sem = nan(numModels,numRGroups);
        for modelIdx = 1:numModels
            R2_indiv(modelIdx,:) = resKept(modelIdx).R2_reg.indiv;
            R2_indiv_sem(modelIdx,:) = resKept(modelIdx).R2_reg_sem.indiv;
        end

        % Determine the number of across-group dimensions required to 
        % predict the other group.
        bestIdxs = nan(1,numRGroups); % Indexes of best models in resKept.
        bestAcross = nan(1,numRGroups); % Across-group dimensionalities of best models
        for groupIdx = 1:numRGroups
            [bestR2, bestIdx] = max(R2_indiv(:,groupIdx));
            best_sem = R2_indiv_sem(bestIdx,groupIdx);
            % NOTE: max used in this manner selects the first instance that
            %       satisfies the input condition.
            [~, bestIdxs(groupIdx)] = max(R2_indiv(:,groupIdx) >= bestR2 - best_sem); 
            bestAcross(groupIdx) = resKept(bestIdxs(groupIdx)).xDim_across;
        end

        % Take the model with the minimum of the optimal values in each direction
        [~, bestGroupIdx] = min(bestAcross);   % Gives index into bestIdxs.
        bestModelIdx = bestIdxs(bestGroupIdx); % Gives index into resKept.
        bestModel = resKept(bestModelIdx);
    end

else
    fprintf('Invalid metric entered.\n')
end
