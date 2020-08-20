function [res, bestModels] = getCrossValResults_fa(runIdx,varargin)
% 
% [res, bestModels] = getCrossValResults_fa(runIdx,...)
%
% Description: Go through cross-validated FA models corresponding to run 
%              runIdx, extract their cross-validated performance 
%              (log-likelihood), and determine the best-performing model.
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
%     res -- (1 x numGroups) cell array; each element is a structure whose 
%            i-th entry has fields
%              xDim   -- int; latent dimensionality
%              sumLL  -- float; cross-validated log-likelihood
%              LL     -- float; average cross-validated log-likelihood
%              LL_sem -- float; standard error of LL across CV folds
%              estParams -- model parameters estimated using all data
%     bestModels -- (1 x numGroups array); indices corresponding to the 
%                   best model for each observation group (based on sumLL)
%                   in res
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 Aug 2020 -- Initial full revision.

baseDir  = '.';
assignopts(who, varargin);
res = [];

runDir = sprintf('%s/mat_results/run%03d', baseDir, runIdx);   
if ~isdir(runDir)
    fprintf('ERROR: %s does not exist.  Exiting...\n', runDir);
    return
else
    D = dir([runDir '/fa*.mat']);
end

if isempty(D)
    fprintf('ERROR: No valid files.  Exiting...\n');
end

for i = 1:length(D)
    P = parseFilename_fa(D(i).name);
    D(i).groupIdx = P.groupIdx;
    D(i).xDim     = P.xDim;
    D(i).cvf      = P.cvf;
end
 
numGroups = max([D.groupIdx]);
% Initialize output structures
res = cell(1,numGroups);
bestModels = nan(1,numGroups);
% Find cross validated performance for each latent dimensionality, and each
% group
for groupIdx = 1:numGroups
    currD = D([D.groupIdx] == groupIdx);
    xDims = unique([currD.xDim]);
    numFolds = max(unique([currD.cvf]));
    numModels = length(xDims);
    for modelIdx = 1:numModels
        xDim = xDims(modelIdx); % Latent dimensionality of current model
        % Initialize output structure
        res{groupIdx}(modelIdx).xDim = xDim;
        % We'll collect performance metrics in the following arrays
        cvLL = nan(numFolds,1);
        fIdxs = find((xDim == [currD.xDim])); % Find files with the appropriate models
        for i = fIdxs    
            fprintf('Loading %s/%s...\n', runDir, currD(i).name);
            ws = load(sprintf('%s/%s', runDir, currD(i).name));
            if ws.cvf == 0
                % For models trained on all data, extract the estimated
                % parameters.
                res{groupIdx}(modelIdx).estParams = ws.estParams;
            else
                % For models trained on CV folds, get performance metrics
                cvLL(ws.cvf) = ws.LLtest;
            end
        end 
        % Compute averages and SEMs
        res{groupIdx}(modelIdx).sumLL = sum(cvLL);
        res{groupIdx}(modelIdx).LL = mean(cvLL);
        res{groupIdx}(modelIdx).LL_sem = std(cvLL,0) / sqrt(numFolds);
    end

    % Find the best model, based on sumLL
    sumLL = [res{groupIdx}.sumLL];
    [~, bestModels(groupIdx)] = max(sumLL);

end