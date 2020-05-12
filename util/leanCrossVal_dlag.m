function res = leanCrossVal_dlag(runIdx, dat, varargin)
%
% res = leanCrossVal_dlag(runIdx, dat, varargin)
% 
% Description: A wrapper function to automate a "lean" approach to
%              selecting DLAG within- and across-group dimensionalities
%              via cross-validation:
%                  (1) Fix all within-group dimensionalities at the 
%                      largest provided value. Sweep over the across-group
%                      dimensionalities provided.
%                  (2) Fix the across-group dimensionality at the value 
%                      given by the best model from step (1). Fix all 
%                      within-group dimensionalities but one at their max
%                      value (if not cross-validated yet) or optimal 
%                      value (if previously cross-validated). Sweep over 
%                      the within-group dimensionalities for the remaining
%                      group. Repeat this step until all optimal 
%                      within-group dimensionalities have been found.
%
% Arguments: 
%
%     NOTE: Arguments and saved results are structured to be the same as 
%           for fit_dlag.m. See fit_dlag.m for details.
%
% Outputs:
%
%     NOTE: Arguments and saved results are structured to be the same as  
%           for fit_dlag.m. See fit_dlag.m for details.
%
%     res -- structure containing saved DLAG results, corresponding to the
%            best-performing model.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     11 May 2020 -- Initial revision.

% Specify defaults for optional arguments
datFormat         = 'seq';
baseDir           = '.';        
overwriteExisting = true; 
saveData          = false;     
method            = 'dlag';      
binWidth          = 20;        
numFolds          = 0;         
yDims             = [];
xDims_across      = [3];
xDims_within      = {[1]};
rGroups           = [1 2];      
startTau          = 2*binWidth;
segLength         = 30;       
init_method       = 'pCCA'; 
learnDelays       = true;   
maxIters          = 1e5;       
freqLL            = 10;          
freqParam         = 100;      
minVarFrac        = 0.01;    
parallelize       = false; 
numWorkers        = 4;
randomSeed        = 0;       
extraOpts         = assignopts(who, varargin);
numGroups         = length(yDims); % Number of groups (areas)

% Fix all within-group dimensionalities at the largest provided value. 
% Sweep over the across-group dimensionalities provided.
curr_xDims_within = cell(1,numGroups);
for groupIdx = 1:numGroups
    curr_xDims_within{groupIdx} = max(xDims_within{groupIdx});
end

fit_dlag(runIdx, dat, ...
         'datFormat', datFormat, ...
         'baseDir', baseDir, ...
         'method', method, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'xDims_across', xDims_across, ...
         'xDims_within', curr_xDims_within, ...
         'yDims', yDims, ...
         'rGroups', rGroups,...
         'startTau', startTau, ...
         'segLength', segLength, ...
         'init_method', init_method, ...
         'learnDelays', learnDelays, ...
         'maxIters', maxIters, ...
         'freqLL', freqLL, ...
         'freqParam', freqParam, ...
         'minVarFrac', minVarFrac, ...
         'parallelize', parallelize, ...
         'randomSeed', randomSeed, ...
         'numWorkers', numWorkers, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);

% Get cross-validation results so far
[cvResults, bestModel] = getCrossValResults_dlag(runIdx, 'baseDir', baseDir);

% Be sure to use the same cross-validation folds from the previous runs
res = getModel_dlag(runIdx, ...
                    cvResults(bestModel).xDim_across, ...
                    cvResults(bestModel).xDim_within, ...
                    'baseDir', baseDir);
if isempty(randomSeed)
    randomSeed = res.rngSettings; 
end

% Fix the across-group dimensionality at the value given by the best model
% from the previous runs. Find the optimal within-group dimensionalities 
% independently.
best_xDim_across = res.xDim_across;
best_xDim_within = nan(1,numGroups); % Collect the optimal within-group dimensionalities 
for groupIdx = 1:numGroups
    curr_xDims_within{groupIdx} = xDims_within{groupIdx};
    fit_dlag(runIdx, dat, ...
             'datFormat', datFormat, ...
             'baseDir', baseDir, ...
             'method', method, ...
             'binWidth', binWidth, ...
             'numFolds', numFolds, ...
             'xDims_across', best_xDim_across, ...
             'xDims_within', curr_xDims_within, ...
             'yDims', yDims, ...
             'rGroups', rGroups,...
             'startTau', startTau, ...
             'segLength', segLength, ...
             'init_method', init_method, ...
             'learnDelays', learnDelays, ...
             'maxIters', maxIters, ...
             'freqLL', freqLL, ...
             'freqParam', freqParam, ...
             'minVarFrac', minVarFrac, ...
             'parallelize', parallelize, ...
             'randomSeed', randomSeed, ...
             'numWorkers', numWorkers, ...
             'overwriteExisting', overwriteExisting, ...
             'saveData', saveData);
         
    % Get cross-validation results so far
    [cvResults, bestModel] = getCrossValResults_dlag(runIdx, 'baseDir', baseDir);
    
    % Fix this within-group dimensionality at its optimal value from here
    % on out.
    curr_xDims_within{groupIdx} = cvResults(bestModel).xDim_within(groupIdx);
    best_xDim_within(groupIdx) = cvResults(bestModel).xDim_within(groupIdx);
    
end

% Return the best-performing model
res = getModel_dlag(runIdx, ...
                    best_xDim_across, ...
                    best_xDim_within, ...
                    'baseDir', baseDir);