function faEngine(seqTrain,seqTest,fname,varargin)
%
% faEngine(seqTrain,seqTest,fname,...)
%
% Description: Initialize and fit FA model parameters; infer latent time
%              courses; and assess generalization performance, if
%              cross-validating. Save results to a file.
%
% Arguments:
%
%     Required:
% 
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
%     fname    -- string; filename where results will be saved
%
%     Optional:
%
%     xDim         -- int; number of across-group dimensions
%     binWidth     -- float; bin width or sample period, in units of time
%                     (default: 20)
%     parallelize  -- logical; Set to true to use Matlab's parfor construct
%                     to parallelize each fold and latent dimensionality 
%                     using multiple cores. (default: false)
%     saveData     -- logical; Set to true to save training and test data
%                     used (could be useful for reproducibility). In
%                     general, NOT RECOMMENDED. Can take up a great deal of
%                     disk space, depending on the dataset used.
%
% Outputs:
%     None. (But saves results to file specified in 'fname')
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

xDim          = 3;
binWidth      = 20;  % in msec
parallelize   = false;
saveData      = false;
extraOpts     = assignopts(who, varargin);
yDim          = size(seqTrain(1).y,1);

% Reformat train data for fitting. Trial ID and location in time do not 
% matter. Static methods treat all data points as independent and 
% identically distributed.
Ytrain = seq2cell2D(seqTrain, yDim, 'datafield', 'y'); % Reusing pCCA code
Ytrain = Ytrain{1};

% =====================
% Fit model parameters
% =====================
if ~parallelize
      fprintf('\nFitting FA model...\n');
end

if xDim == 0
    % If xDim = 0, that means we're checking to see if there's
    % no interaction between observed variables. 
    estParams = indepGaussFit(Ytrain);
    [LLtrain, ~] = indepGaussEval(Ytrain, estParams);
else
    % Otherwise, fit FA normally.
    [estParams, LL] = em_fa(Ytrain, xDim, ...
                            'parallelize', parallelize, extraOpts{:});
    
    % Extract neural trajectories using learned parameters
    [Xtrain, LLtrain] = fa_estep(Ytrain, estParams);
    
    % Reformat inferred latents according to trials and add to seqTrain
    seqTrain = segmentByTrial(seqTrain, Xtrain.mean, 'xsm');
end

% ==================================
% Assess generalization performance
% ==================================
if ~isempty(seqTest)
    % Reformat test data, as was done for train data
    Ytest = seq2cell2D(seqTest, yDim, 'datafield', 'y'); % Reusing pCCA code
    Ytest = Ytest{1};
    
    if xDim == 0
        [LLtest, ~] = indepGaussEval(Ytest, estParams);
    
    else
        % Log-likelihood of test data
        [Xtest, LLtest] = fa_estep(Ytest, estParams);

        % Reformat inferred latents according to trials and add to seqTest
        seqTest = segmentByTrial(seqTest, Xtest.mean, 'xsm');
    end
end

% ===========
% Save results
% ===========

vars = who;
if ~parallelize
      fprintf('Saving %s...\n', fname);
end

% Remove variables we don't want before saving.
vars = vars(~ismember(vars,{'Ytrain', 'Xtrain', 'Ytest', 'Xtest'}));
% Remove redundant variables before saving.
vars = vars(~ismember(vars, {'varargin'}));
% If saveData is true, then we'll write seqTrain, seqTest. 
% Otherwise, we won't.
if ~saveData
    vars = vars(~ismember(vars, {'seqTrain', 'seqTest'}));
end
save(fname,vars{:});
