function [res, bestModel] = getCrossValResults_dlag(runIdx,varargin)
% 
% [res, bestModel] = getCrossValResults_dlag(runIdx,...)
%
% Description: Go through cross-validated DLAG models corresponding to run 
%              runIdx, extract their cross-validated performance, and
%              determine the best-performing model.
% Arguments:
%
%     Required:
%
%     runIdx    -- int; results files will be loaded from 
%                  baseDir/mat_results/runXXX, where XXX is runIdx.
%                  baseDir can be specified by the user (see below)
% 
%     Optional:
%
%     baseDir   -- string; specifies directory in which to store
%                  mat_results. (default: '.', i.e., current directory)
%     verbose   -- logical; set true to print out which files are being
%                  accessed (default: true)
% 
% Outputs:
%
%     res -- structure whose i-th entry has fields
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
%              R2orth_reg.indiv -- (xDim_across+1 x 2) array; average 
%                                  cross-validated R^2 error in each 
%                                  direction for orthonormalized DLAG 
%                                  predictions, for the pair in rGroups.
%              R2orth_reg_sem.indiv -- (xDim_across+1 x 2) array; standard error 
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
%              MSEorth_reg.indiv -- (xDim_across+1 x 2) array; average 
%                                   cross-validated mean-squared error in 
%                                   each direction for orthonormalized DLAG 
%                                   predictions, for the pair in rGroups.
%              MSEorth_reg_sem.indiv -- (xDim_across+1 x 2) array; standard 
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
%     bestModel -- int; index corresponding to the best model 
%                  (based on sumLL) in res
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     09 Apr 2020 -- Initial full revision.
%     23 May 2020 -- Added warning for models with errors during fitting.
%     11 Jun 2020 -- Updated for expanded cross-validation metrics.
%     12 Mar 2022 -- Added verbose option.

baseDir  = '.';
verbose = true;
assignopts(who, varargin);
res = [];
bestModel = [];

runDir = sprintf('%s/mat_results/run%03d', baseDir, runIdx);   
if ~isdir(runDir)
    fprintf('ERROR: %s does not exist.  Exiting...\n', runDir);
    return
else
    D = dir([runDir '/dlag*.mat']);
end

if isempty(D)
    fprintf('ERROR: No valid files.  Exiting...\n');
    return;
end

for i = 1:length(D)
    P = parseFilename_dlag(D(i).name);
    D(i).xDim_across = P.xDim_across;
    D(i).xDim_within = P.xDim_within;
    D(i).xDim_all = {[D(i).xDim_across; D(i).xDim_within']};
    D(i).cvf = P.cvf;
end

% Find cross validated performance for each combination of within- and
% across-group latent dimensionalities
xDims_all = unique(cell2mat([D.xDim_all])', 'rows');
numFolds = max(unique([D.cvf]));
numModels = size(xDims_all,1);
clear res;
for modelIdx = 1:numModels
    xDim_across = xDims_all(modelIdx,1); % Across-group dimensionality of current model
    xDim_within = xDims_all(modelIdx,2:end); % Within-group dimensionalities
    xDim_all = xDim_across + xDim_within;
    numGroups = length(xDim_within);
    
    % Initialize output structure
    res(modelIdx).xDim_across = xDim_across;
    res(modelIdx).xDim_within = xDim_within;
    
    % We'll collect performance metrics in the following arrays
    % Log-likelihood, joint and individual
    cvLL.joint = nan(numFolds,1);
    cvLL.indiv = nan(numFolds,numGroups);
    
    % R^2 and MSE for pairwise regression, joint and individual
    cvR2_reg.joint = nan(numFolds,1);
    cvR2_reg.indiv = nan(numFolds,2);
    cvR2orth_reg.indiv = nan(numFolds,xDim_across+1,2);
    
    cvMSE_reg.joint = nan(numFolds,1);
    cvMSE_reg.indiv = nan(numFolds,2);
    cvMSEorth_reg.indiv = nan(numFolds,xDim_across+1,2);
    
    % R^2 and MSE for denoised reconstruction, joint and individual
    cvR2_denoise.joint = nan(numFolds,1);
    cvMSE_denoise.joint = nan(numFolds,1);
    
    cvR2_denoise.indiv = nan(numFolds,numGroups);
    cvMSE_denoise.indiv = nan(numFolds,numGroups);
    
    % Orthonormalized DLAG model performance
    cvR2orth_denoise.indiv = cell(1,numGroups);
    cvMSEorth_denoise.indiv = cell(1,numGroups);
    for groupIdx = 1:numGroups
        cvR2orth_denoise.indiv{groupIdx} = nan(numFolds,xDim_all(groupIdx)+1);
        cvMSEorth_denoise.indiv{groupIdx} = nan(numFolds,xDim_all(groupIdx)+1);
    end
    
    % Find files with the appropriate models
    fIdxs = find(ismember(cell2mat([D.xDim_all])', xDims_all(modelIdx,:), 'rows'))'; 
    for i = fIdxs 
        if verbose
            fprintf('Loading %s/%s...\n', runDir, D(i).name);
        end
        ws = load(sprintf('%s/%s', runDir, D(i).name));
        if ws.err_status
            % Flag models where fitting stopped due to an error
            fprintf('Warning: Data likelihood decreased for %s\n', D(i).name);
        end
        if ws.cvf == 0
            % For models trained on all data, extract the estimated
            % parameters.
            res(modelIdx).estParams = ws.estParams;
        else
            % For models trained on CV folds, get performance metrics
            
            % Log-likelihood, joint and individual
            cvLL.joint(ws.cvf) = ws.LLtest.joint;
            cvLL.indiv(ws.cvf,:) = ws.LLtest.indiv;

            % R^2 and MSE for pairwise regression, joint and individual
            cvR2_reg.joint(ws.cvf) = ws.R2_reg.joint;
            cvR2_reg.indiv(ws.cvf,:) = ws.R2_reg.indiv;
            cvR2orth_reg.indiv(ws.cvf,:,:) = ws.R2orth_reg.indiv;

            cvMSE_reg.joint(ws.cvf) = ws.MSE_reg.joint;
            cvMSE_reg.indiv(ws.cvf,:) = ws.MSE_reg.indiv;
            cvMSEorth_reg.indiv(ws.cvf,:,:) = ws.MSEorth_reg.indiv;

            % R^2 and MSE for denoised reconstruction, joint and individual
            cvR2_denoise.joint(ws.cvf) = ws.R2_denoise.joint;
            cvMSE_denoise.joint(ws.cvf) = ws.MSE_denoise.joint;

            cvR2_denoise.indiv(ws.cvf,:) = [ws.R2_denoise.indiv];
            cvMSE_denoise.indiv(ws.cvf,:) = [ws.MSE_denoise.indiv];

            % Orthonormalized DLAG model performance
            for groupIdx = 1:numGroups
                cvR2orth_denoise.indiv{groupIdx}(ws.cvf,:) = ws.R2orth_denoise.indiv{groupIdx};
                cvMSEorth_denoise.indiv{groupIdx}(ws.cvf,:) = ws.MSEorth_denoise.indiv{groupIdx};
            end
        end
    end 
    % Compute averages and SEMs
    res(modelIdx).rGroups = ws.rGroups;

    % Log-likelihood, joint and individual
    res(modelIdx).sumLL.joint = sum(cvLL.joint);
    res(modelIdx).LL.joint = mean(cvLL.joint);
    res(modelIdx).LL_sem.joint = std(cvLL.joint,0) / sqrt(numFolds);
    res(modelIdx).sumLL.indiv = sum(cvLL.indiv,1);
    res(modelIdx).LL.indiv = mean(cvLL.indiv,1);
    res(modelIdx).LL_sem.indiv = std(cvLL.indiv,0,1) /sqrt(numFolds);
    
    % R^2 and MSE for pairwise regression, joint and individual
    res(modelIdx).R2_reg.joint = mean(cvR2_reg.joint);
    res(modelIdx).R2_reg_sem.joint = std(cvR2_reg.joint,0) ./ sqrt(numFolds);
    res(modelIdx).R2_reg.indiv = mean(cvR2_reg.indiv,1);
    res(modelIdx).R2_reg_sem.indiv = std(cvR2_reg.indiv,0,1) ./ sqrt(numFolds);
    res(modelIdx).R2orth_reg.indiv = reshape(mean(cvR2orth_reg.indiv,1), [xDim_across+1,2,1]);
    res(modelIdx).R2orth_reg_sem.indiv = reshape(std(cvR2orth_reg.indiv,0,1) ./ sqrt(numFolds), [xDim_across+1,2,1]);

    res(modelIdx).MSE_reg.joint = mean(cvMSE_reg.joint);
    res(modelIdx).MSE_reg_sem.joint = std(cvMSE_reg.joint,0) ./ sqrt(numFolds);
    res(modelIdx).MSE_reg.indiv = mean(cvMSE_reg.indiv,1);
    res(modelIdx).MSE_reg_sem.indiv = std(cvMSE_reg.indiv,0,1) ./ sqrt(numFolds);
    res(modelIdx).MSEorth_reg.indiv = reshape(mean(cvMSEorth_reg.indiv,1), [xDim_across+1,2,1]);
    res(modelIdx).MSEorth_reg_sem.indiv = reshape(std(cvMSEorth_reg.indiv,0,1) ./ sqrt(numFolds), [xDim_across+1,2,1]);

    % R^2 and MSE for denoised reconstruction, joint and individual
    res(modelIdx).R2_denoise.joint = mean(cvR2_denoise.joint);
    res(modelIdx).R2_denoise_sem.joint = std(cvR2_denoise.joint,0) ./ sqrt(numFolds);

    res(modelIdx).MSE_denoise.joint = mean(cvMSE_denoise.joint);
    res(modelIdx).MSE_denoise_sem.joint = std(cvMSE_denoise.joint,0) ./ sqrt(numFolds);

    res(modelIdx).R2_denoise.indiv = mean(cvR2_denoise.indiv,1);
    res(modelIdx).R2_denoise_sem.indiv = std(cvR2_denoise.indiv,0,1) ./ sqrt(numFolds);

    res(modelIdx).MSE_denoise.indiv = mean(cvMSE_denoise.indiv,1);
    res(modelIdx).MSE_denoise_sem.indiv = std(cvMSE_denoise.indiv,0,1) ./ sqrt(numFolds);
    
    % Orthonormalized DLAG model performance
    for groupIdx = 1:numGroups
        res(modelIdx).R2orth_denoise.indiv{groupIdx} = mean(cvR2orth_denoise.indiv{groupIdx},1);
        res(modelIdx).R2orth_denoise_sem.indiv{groupIdx} = std(cvR2orth_denoise.indiv{groupIdx},0,1) ./ sqrt(numFolds);

        res(modelIdx).MSEorth_denoise.indiv{groupIdx} = mean(cvMSEorth_denoise.indiv{groupIdx},1);
        res(modelIdx).MSEorth_denoise_sem.indiv{groupIdx} = std(cvMSEorth_denoise.indiv{groupIdx},0,1) ./ sqrt(numFolds);
    end
end

% Find the best model, based on sumLL.joint
sumLL = nan(1,numModels);
for modelIdx = 1:numModels
    sumLL(modelIdx) = res(modelIdx).sumLL.joint;
end
[~, bestModel] = max(sumLL);