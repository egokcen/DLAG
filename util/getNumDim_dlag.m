function dim = getNumDim_dlag(res,varargin)
% 
% dim = getNumDim_dlag(res,...)
%
% Description: Determine the number of across-group dimensions, and the
%              number of overall dimensions required to explain each area,
%              based on the cross-validated results in res.
%
% Arguments:
%
%     Required:
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
%     Optional:
%
%     metric -- string; choose which cross-validated performance metric to
%               use for selection ('MSE', 'R2') (default: 'R2')
% 
% Outputs:
%     dim -- structure with the following fields
%            all    -- (1 x numGroups) array; Overall number of dimensions 
%                      (within and across) required to explain each group
%            across -- int; Number of across-group dimensions required
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Jun 2020 -- Initial full revision.

metric = 'R2';
assignopts(who, varargin);

numGroups = length(res.xDim_within);
rGroups = res.rGroups;

% Initialize output structure
dim.all = nan(1,numGroups);
dim.across = nan(1,numGroups);

if isequal(metric,'MSE')
    for groupIdx = 1:numGroups
        % Determine the number of overall dimensions required to explain 
        % the current group
        MSE_denoise = res.MSEorth_denoise.indiv{groupIdx};
        MSE_denoise_sem = res.MSEorth_denoise_sem.indiv{groupIdx};
        [bestMSE, bestIdx] = min(MSE_denoise);
        best_sem = MSE_denoise_sem(bestIdx);
        % NOTE: max used in this manner selects the first instance that
        %       satisfies the input condition.
        [~, dim.all(groupIdx)] = max(MSE_denoise <= bestMSE + best_sem);
        
        % Determine the number of across-group dimensions required to 
        % predict the other group.
        if ismember(groupIdx,rGroups)
            sourceIdx = find(rGroups == groupIdx);
            MSE_reg = res.MSEorth_reg.indiv(:,sourceIdx);
            MSE_reg_sem = res.MSEorth_reg_sem.indiv(:,sourceIdx);
            [bestMSE, bestIdx] = min(MSE_reg);
            best_sem = MSE_reg_sem(bestIdx);
            % NOTE: max used in this manner selects the first instance that
            %       satisfies the input condition.
            [~, dim.across(groupIdx)] = max(MSE_reg <= bestMSE + best_sem);
        end
    end
    
    % Determine the final number of across-group dimensions required
    dim.across = min(dim.across);

elseif isequal(metric,'R2')
    for groupIdx = 1:numGroups
        % Determine the number of overall dimensions required to explain 
        % the current group
        R2_denoise = res.R2orth_denoise.indiv{groupIdx};
        R2_denoise_sem = res.R2orth_denoise_sem.indiv{groupIdx};
        [bestR2, bestIdx] = max(R2_denoise);
        best_sem = R2_denoise_sem(bestIdx);
        % NOTE: max used in this manner selects the first instance that
        %       satisfies the input condition.
        [~, dim.all(groupIdx)] = max(R2_denoise >= bestR2 - best_sem);
        
        % Determine the number of across-group dimensions required to 
        % predict the other group.
        if ismember(groupIdx,rGroups)
            sourceIdx = find(rGroups == groupIdx);
            R2_reg = res.R2orth_reg.indiv(:,sourceIdx);
            R2_reg_sem = res.R2orth_reg_sem.indiv(:,sourceIdx);
            [bestR2, bestIdx] = max(R2_reg);
            best_sem = R2_reg_sem(bestIdx);
            % NOTE: max used in this manner selects the first instance that
            %       satisfies the input condition.
            [~, dim.across(groupIdx)] = max(R2_reg >= bestR2 - best_sem);
        end
    end
    
    % Determine the final number of across-group dimensions required
    dim.across = min(dim.across);

else
    fprintf('Invalid metric entered.\n')
end
