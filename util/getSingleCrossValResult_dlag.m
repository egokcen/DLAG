function res = getSingleCrossValResult_dlag(cvResults,xDim_across,xDim_within)
% 
% res = getSingleCrossValResult_dlag(cvResults,xDim_across,xDim_within)
%
% Description: Access the cross-validation results for a single model in
%              the cvResults structure.
% Arguments:
%
%     cvResults -- structure whose i-th entry has fields
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
%              R2orth_reg.indiv -- (xDim_across x 2) array; average 
%                                  cross-validated R^2 error in each 
%                                  direction for orthonormalized DLAG 
%                                  predictions, for the pair in rGroups.
%              R2orth_reg_sem.indiv -- (xDim_across x 2) array; standard error 
%                                      of R2orth_reg.indiv across CV folds
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
%              MSEorth_reg.indiv -- (xDim_across x 2) array; average 
%                                   cross-validated mean-squared error in 
%                                   each direction for orthonormalized DLAG 
%                                   predictions, for the pair in rGroups.
%              MSEorth_reg_sem.indiv -- (xDim_across x 2) array; standard 
%                                       error of MSEorth_reg.indiv across
%                                       CV folds
%              R2_denoise.joint -- float; average cross-validated R^2,
%                                  evaluated jointly from denoised 
%                                  reconstruction across all groups
%              R2_denoise.indiv -- (1 x numGroups) array; average 
%                                  cross-validated R^2, evaluated for each 
%                                  group from denoised reconstruction
%              R2_denoise_sem.joint -- float; standard error of
%                                      R2_denoise.joint across CV folds
%              R2_denoise_sem.indiv -- (1 x numGroups) array; standard error 
%                                      of R2_denoise.indiv across CV folds
%              R2orth_denoise.indiv -- (1 x numGroups) cell array;
%                                      average cross-validated R^2 for the 
%                                      orthonormalized DLAG model, evaluated
%                                      for each group from denoised
%                                      reconstruction
%              R2orth_denoise_sem.indiv -- (1 x numGroups) cell array; 
%                                      standard error of 
%                                      R2orth_denoise.indiv across CV folds
%              MSE_denoise.joint -- float; average cross-validated 
%                                   mean-squared error, evaluated jointly
%                                   from denoised reconstruction across 
%                                   all groups
%              MSE_denoise.indiv -- (1 x numGroups) array; average 
%                                   cross-validated mean-squared error, 
%                                   evaluated for each group from denoised
%                                   reconstruction
%              MSE_denoise_sem.joint -- float; standard error of
%                                       MSE_denoise.joint across CV folds
%              MSE_denoise_sem.indiv -- (1 x numGroups) array; standard 
%                                       error of MSE_denoise.indiv across 
%                                       CV folds
%              MSEorth_denoise.indiv -- (1 x numGroups) cell array;
%                                       average cross-validated mean-squared 
%                                       error for the orthonormalized DLAG 
%                                       model, evaluated for each group 
%                                       from denoised reconstruction
%              MSEorth_denoise_sem.indiv -- (1 x numGroups) cell array; 
%                                      standard error of 
%                                      MSEorth_denoise.indiv across CV folds
%              estParams -- model parameters estimated using all data
%              rGroups -- (1 x 2) array; the indexes of two groups 
%                         used to measure generalization performance
%                         via regression
%
%     xDim_across -- int; Number of across-group dimensions
%     xDim_within -- (1 x numGroups) array; Number of within-group
%                    dimensions in each group
% 
% Outputs:
%
%     res -- struct; The entry of cvResults corresponding the model with 
%            specified within- and across-group dimensions.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     20 Jun 2020 -- Initial full revision.

res = [];
xDims_desired = [xDim_across xDim_within];
xDims_list = [vertcat(cvResults.xDim_across) vertcat(cvResults.xDim_within)];
[modelPresent, modelIdx] = ismember(xDims_desired, xDims_list, 'rows');
if ~modelPresent
    fprintf('ERROR: Requested model not present in given results structure.\n');
else
    res = cvResults(modelIdx);
end
