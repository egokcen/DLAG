function [res, bestModel] = getCrossValResults_pcca(runIdx,varargin)
% 
% [res, bestModel] = getCrossValResults_pcca(runIdx,...)
%
% Description: Go through cross-validated pCCA models corresponding to run 
%              runIdx, extract their cross-validated performance 
%              (log-likelihood; MSE and R^2 for pairwise regression), and
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
% 
% Outputs:
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
%     bestModel -- int; index corresponding to the best model 
%                  (based on sumLL) in res
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     09 Apr 2020 -- Initial full revision.

baseDir  = '.';
assignopts(who, varargin);
res = [];

runDir = sprintf('%s/mat_results/run%03d', baseDir, runIdx);   
if ~isdir(runDir)
    fprintf('ERROR: %s does not exist.  Exiting...\n', runDir);
    return
else
    D = dir([runDir '/pcca*.mat']);
end

if isempty(D)
    fprintf('ERROR: No valid files.  Exiting...\n');
end

for i = 1:length(D)
    P = parseFilename_pcca(D(i).name);
    D(i).xDim   = P.xDim;
    D(i).cvf    = P.cvf;
end
  
% Find cross validated performance for each latent dimensionality
xDims = unique([D.xDim]);
numFolds = max(unique([D.cvf]));
numModels = length(xDims);
clear res;
for modelIdx = 1:numModels
    xDim = xDims(modelIdx); % Latent dimensionality of current model
    % Initialize output structure
    res(modelIdx).xDim = xDim;
    % We'll collect performance metrics in the following arrays
    cvLL = nan(numFolds,1);
    cvR2 = nan(numFolds,2);
    cvMSE = nan(numFolds,2);
    fIdxs = find((xDim == [D.xDim])); % Find files with the appropriate models
    for i = fIdxs    
        fprintf('Loading %s/%s...\n', runDir, D(i).name);
        ws = load(sprintf('%s/%s', runDir, D(i).name));
        if ws.cvf == 0
            % For models trained on all data, extract the estimated
            % parameters.
            res(modelIdx).estParams = ws.estParams;
            res(modelIdx).rGroups = ws.rGroups;
        else
            % For models trained on CV folds, get performance metrics
            cvLL(ws.cvf) = ws.LLtest;
            cvR2(ws.cvf,:) = ws.R2;
            cvMSE(ws.cvf,:) = ws.MSE;
        end
    end 
    % Compute averages and SEMs
    res(modelIdx).sumLL = sum(cvLL);
    res(modelIdx).LL = mean(cvLL);
    res(modelIdx).LL_sem = std(cvLL,0) / sqrt(numFolds);
    res(modelIdx).R2 = mean(cvR2,1);
    res(modelIdx).R2_sem = std(cvR2,0,1) / sqrt(numFolds);
    res(modelIdx).MSE = mean(cvMSE,1);
    res(modelIdx).MSE_sem = std(cvMSE,0,1) / sqrt(numFolds);
end

% Find the best model, based on sumLL
sumLL = [res.sumLL];
[~, bestModel] = max(sumLL);

