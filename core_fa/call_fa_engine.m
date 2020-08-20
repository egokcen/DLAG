function call_fa_engine(fname,seqTrain,seqTest,varargin)
%
% call_fa_engine(fname,seqTrain,seqTest,...)
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
%     method   -- string; method to be used (right now, just 'fa') 
%                 (default: 'fa')
%     binWidth -- float; bin width or sample period, in unitsof time
%                 (default: 20)
%     cvf      -- int; indicates which cross-validation fold input data
%                 corresponds to(default: 0).
%                 0 indicates no cross-validation, i.e. train on all trials.
%     xDim     -- int; number of latent dimensions
%     groupIdx -- int; indicates which observation group is being fitted.
%                 (default: 1)
%     overwriteExisting -- logical; Set to true to overwrite existing
%                          results files, if their names match the the ones
%                          created during this run. (default: true)
%     rngSettings -- structure with the random number generator 
%                    settings used during run time. Includes 
%                    fields 'Type', 'Seed', and 'State'.
%
% Outputs:
%
%     Results are saved to a file. See fit_fa.m for details.
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

method        = 'fa';
binWidth      = 20; % in msec
cvf           = 0;
xDim          = 3;
groupIdx      = 1;
overwriteExisting = true;
rngSettings   = [];
extraOpts     = assignopts(who, varargin);

extraOpts = {extraOpts{:},'cvf',cvf}; % cvf will continue to get passed on

% Skip existing files, if overwriting is not desired
if ~(overwriteExisting) && exist([fname '.mat'], 'file')
    fprintf('%s already exists.  Skipping...\n', fname);
    return;
end

% The following does the heavy lifting for learning and inference.
if isequal(method,'fa')
    faEngine(seqTrain, seqTest, fname, ...
        'xDim', xDim, 'binWidth', binWidth, extraOpts{:});
end

% Save results to a file.
if exist([fname '.mat'], 'file')
    save(fname, 'method', 'cvf', 'rngSettings', 'groupIdx', '-append');
end
