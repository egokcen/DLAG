function result = call_dlag_engine(fname,seqTrain,seqTest,varargin)
%
% result = call_dlag_engine(fname,seqTrain,seqTest,...)
%
% Description: Interfaces with DLAG Engines, which do the heavy lifting for
%              learning and inference.
%              NOTE: For now, only dlagEngine is supported in this
%                    function, but it is structured to handle multiple new
%                    versions in the future.
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
%     method   -- string; method to be used (currently one supported):
%                 'dlag'. (default: 'dlag')
%     binWidth -- float; bin width or sample period, in units of time
%                 (default: 20)
%     cvf      -- int; indicates which cross-validation fold input data
%                 corresponds to(default: 0).
%                 0 indicates no cross-validation, i.e. train on all trials.
%     xDim_across -- int; number of across-group dimensions
%     xDim_within -- (1 x numGroups) array; number of within-
%                    group dimensions for each group
%     overwriteExisting -- logical; Set to true to overwrite existing
%                          results files, if their names match the the ones
%                          created during this run. (default: true)
%     rngSettings -- structure with the random number generator 
%                    settings used during run time. Includes 
%                    fields 'Type', 'Seed', and 'State'.
%
% Outputs:
%
%     result -- structure containing all variables saved in 
%               ./mat_results/runXXX/ if 'numFolds' is 0.  
%               Else, the structure is empty. See fit_dlag.m for details.
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 Mar 2020 -- Initial full revision.
%     07 May 2020 -- Save random number generator settings.
%     13 May 2020 -- Removed hasSpikesBool functionality.
%     24 May 2020 -- Printed info about fitted model moved to fit_dlag.m

method        = 'dlag';
binWidth      = 20; % in msec
cvf           = 0;
xDim_across   = 3;
xDim_within   = [];
overwriteExisting = true;
rngSettings   = [];
extraOpts     = assignopts(who, varargin);
numGroups     = length(xDim_within); % Number of groups (areas)

% If doing cross-validation, don't use private noise variance floor.
if cvf > 0
    extraOpts = {extraOpts{:}, 'minVarFrac', -Inf};
end
extraOpts = {extraOpts{:},'cvf',cvf}; % cvf will continue to get passed on

% Skip existing files, if overwriting is not desired
if ~(overwriteExisting) && exist([fname '.mat'], 'file')
    fprintf('%s already exists.  Skipping...\n', fname);
    return;
end

% The following does the heavy lifting for learning and inference.
if isequal(method,'dlag')
    dlagEngine(seqTrain, seqTest, fname,...
        'xDim_across', xDim_across, 'xDim_within', xDim_within, ...
        'binWidth', binWidth, extraOpts{:});
end

if exist([fname '.mat'], 'file')
    save(fname, 'method', 'cvf', 'rngSettings', '-append');
end

% Results are saved to a file.
% Here, return those results (from the file) only if not doing
% cross-validation
result = {};
if (nargout == 1) && (numFolds == 0) && exist([fname '.mat'], 'file')
    result = load(fname);
end