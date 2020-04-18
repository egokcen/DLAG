function result = call_pcca_engine(fname,seqTrain,seqTest,varargin)
%
% result = call_pcca_engine(fname,seqTrain,seqTest,...)
%
% Description: Interfaces with static method engines, which do the heavy
%              lifting for learning and inference.
%
% Arguments:
%
%     Required:
% 
%     fname    -- string; filename where results will be saved
%     seqTrain -- training data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
% 
%     seqTest  -- testing data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
%     Optional:
%
%     method   -- string; method to be used (right now, just 'pcca') 
%                 (default: 'pcca)
%     binWidth -- float; bin width or sample period, in unitsof time
%                 (default: 20)
%     cvf      -- int; indicates which cross-validation fold input data
%                 corresponds to(default: 0).
%                 0 indicates no cross-validation, i.e. train on all trials.
%     xDim     -- int; number of across-group dimensions
%     yDims    -- (1 x numGroups) array; Specify the number features 
%                 (neurons) in each group (area). Elements in yDims
%                 should match the format of data in seqTrain and seqTest. 
%     hasSpikesBool -- (1 x yDim) logical array; Each entry is 1 if an
%                      observation (neuron) had any non-zero activity, 0
%                      otherwise. (default: [])
%     overwriteExisting -- logical; Set to true to overwrite existing
%                          results files, if their names match the the ones
%                          created during this run. (default: true)
%
% Outputs:
%
%     result -- structure containing all variables saved in 
%               ./mat_results/runXXX/ if 'numFolds' is 0.  
%               Else, the structure is empty. See fit_pcca.m for details.
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 Mar 2020 -- Initial full revision.

method        = 'pcca';
binWidth      = 20; % in msec
cvf           = 0;
yDims         = [];
xDim          = 3;
hasSpikesBool = [];
overwriteExisting = true;
extraOpts     = assignopts(who, varargin);

extraOpts = {extraOpts{:},'cvf',cvf}; % cvf will continue to get passed on

% Skip existing files, if overwriting is not desired
if ~(overwriteExisting) && exist([fname '.mat'], 'file')
    fprintf('%s already exists.  Skipping...\n', fname);
    return;
end

% Print useful info about the model about to be fitted.
fprintf('Number of training trials: %d\n', length(seqTrain));
fprintf('Number of test trials: %d\n', length(seqTest));
fprintf('Across-group latent dimensionality: %d\n', xDim);
fprintf('Observation dimensionality: %d\n', size(seqTrain(1).y,1));

% The following does the heavy lifting for learning and inference.
if isequal(method,'pcca')
    pccaEngine(seqTrain, seqTest, fname, 'yDims', yDims, ...
        'xDim', xDim, 'binWidth', binWidth, extraOpts{:});
end

if exist([fname '.mat'], 'file')
    save(fname, 'method', 'cvf', 'hasSpikesBool', '-append');
end

% Results are saved to a file.
% Here, return those results (from the file) only if not doing
% cross-validation
result = {};
if (nargout == 1) && (numFolds == 0) && exist([fname '.mat'], 'file')
    result = load(fname);
end