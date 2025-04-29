function pccaEngine(seqTrain,seqTest,fname,varargin)
%
% pccaEngine(seqTrain,seqTest,fname,...)
%
% Description: Initialize and fit pCCA model parameters; infer latent time
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
%     yDims        -- (1 x numGroups) array; Specify the number features 
%                     (neurons) in each group (area). Elements in yDims
%                     should match the format of data in seqTrain and seqTest. 
%     binWidth     -- float; bin width or sample period, in units of time
%                     (default: 20)
%     rGroups      -- (1 x 2) array; Used to assess cross-validated
%                     performance via pairwise regression. Each element
%                     specifies a group to be included in the regression.
%                     (default: [1 2])
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
%     18 Mar 2020 -- Initial full revision.

xDim          = 3;
yDims         = [];
rGroups       = [1 2];
binWidth      = 20;  % in msec
parallelize   = false;
saveData      = false;
extraOpts     = assignopts(who, varargin);

% Reformat train data for fitting. Trial ID and location in time do not 
% matter. Static methods treat all data points as independent and 
% identically distributed.
Ys_train = seq2cell2D(seqTrain, yDims, 'datafield', 'y');

% =====================
% Fit model parameters
% =====================
if ~parallelize
      fprintf('\nFitting pCCA model...\n');
end

if xDim == 0
    % If xDim = 0, that means we're checking to see if there's
    % no interaction between the groups. 
    estParams = indepGroupFit(Ys_train);
    [LLtrain, ~, ~] = indepGroupEval(Ys_train, estParams, rGroups);
else
    % Otherwise, fit pCCA normally.
    [estParams, LL] = em_pcca(Ys_train, xDim, ...
                              'parallelize', parallelize, extraOpts{:});
    
    % Extract neural trajectories using learned parameters
    [Xtrain, LLtrain] = pcca_estep(Ys_train, estParams);
    
    % Reformat inferred latents according to trials and add to seqTrain
    seqTrain = segmentByTrial(seqTrain, Xtrain.mean, 'xsm');
end

% ==================================
% Assess generalization performance
% ==================================
if ~isempty(seqTest)
    % Reformat test data, as was done for train data
    Ys_test = seq2cell2D(seqTest, yDims, 'datafield', 'y');
    
    if xDim == 0
        [LLtest, R2, MSE] = indepGroupEval(Ys_test, estParams, rGroups);
    
    else
        % Pairwise regression on test data
        [~, R2, MSE] = pairwise_regress_pcca(Ys_test, estParams, rGroups);

        % Log-likelihood of test data
        [Xtest, LLtest] = pcca_estep(Ys_test, estParams);

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
vars = vars(~ismember(vars,{'Ys_train', 'Xtrain', 'Ys_test', 'Xtest'}));
% Remove redundant variables before saving.
vars = vars(~ismember(vars, {'varargin'}));
% If saveData is true, then we'll write seqTrain, seqTest. 
% Otherwise, we won't.
if ~saveData
    vars = vars(~ismember(vars, {'seqTrain', 'seqTest'}));
end
save(fname,vars{:});
