function dlagEngine(seqTrain,seqTest,fname,varargin)
%
% dlagEngine(seqTrain,seqTest,fname,...)
%
% Description: Initialize and fit DLAG model parameters; infer latent time
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
%     xDim_across  -- int; number of across-group dimensions
%     xDim_within  -- (1 x numGroups) array; number of within-
%                     group dimensions for each group
%     yDims        -- (1 x numGroups) array; Specify the number features 
%                     (neurons) in each group (area). Elements in yDims
%                     should match the format of data in seqTrain and seqTest. 
%     binWidth     -- float; bin width or sample period, in units of time
%                     (default: 20)
%     startTau     -- float; Initial GP timescale, in units of time (same
%                     binWidth) (default: 2*binWidth)
%     startEps     -- float; Initial GP noise variance (default: 1e-3)
%     startDelay   -- (1 x numGroups) array; Initial delays between 
%                     across-area latents and groups. All across-area 
%                     latents reach a particular group after the same
%                     specified delay. Different groups can have different
%                     delays. Entries in units of time (same as binWidth). 
%                     (default: [])
%     rGroups      -- (1 x 2) array; Used to assess cross-validated
%                     performance via pairwise regression. Each element
%                     specifies a group to be included in the regression.
%                     (default: [1 2])
%     init_method  -- string; Specify how DLAG parameters should be
%                     initialized. Options currently supported:
%                         'pCCA'   -- parameters are initialized using 
%                                     probabilistic canonical correlation 
%                                     analysis (pCCA), described in 
%                                     (Bach & Jordan, 2005).
%                         'params' -- parameters are initialized to a 
%                                     particular realization of DLAG
%                                     model parameters. For example, one 
%                                     may wish to continue training
%                                     a DLAG model for which training was 
%                                     cut short, or performed in a 
%                                     different context.
%                     (default: 'pCCA')
%     covType      -- string; Specify GP covariance kernel type. Options
%                     currently supported:
%                         'rbf' -- Radial basis function, or squared
%                                  exponential kernel
%                     (default: 'rbf')
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
%     09 Jun 2020 -- Expanded metrics of generalization performance.

xDim_across   = 3;
xDim_within   = [];
yDims         = [];
rGroups       = [1 2];
binWidth      = 20;  
startTau      = 2*binWidth; 
startEps      = 1e-3;
startDelay    = [];  
init_method   = 'pCCA';
covType       = 'rbf';
parallelize   = false;
saveData      = false;
extraOpts     = assignopts(who, varargin);
numGroups     = length(yDims);

% For compute efficiency, train on equal-length segments of trials
seqTrainCut = cutTrials(seqTrain, extraOpts{:});
if isempty(seqTrainCut)
    fprintf('WARNING: no segments extracted for training.  Defaulting to segLength=Inf.\n');
    seqTrainCut = cutTrials(seqTrain, 'segLength', Inf);
end

% ========================================
% Initialize model parameters
% ========================================
startParams = initialize_dlag(seqTrainCut, init_method, ...
    'binWidth', binWidth, 'covType', covType, 'xDim_across', xDim_across, ...
    'xDim_within', xDim_within, 'startEps', startEps, 'startTau', startTau, ...
    'yDims', yDims, 'rGroups', rGroups, 'startDelay', startDelay, ...
    'parallelize', parallelize, extraOpts{:});

% =====================
% Fit model parameters
% =====================
if ~parallelize
      fprintf('\nFitting DLAG model...\n');
end

[estParams, seqTrainCut, LLcut, iterTime, D, gams_across, gams_within, err_status, msg] ...
    = em_dlag(startParams, seqTrainCut, ...
                  'parallelize', parallelize, extraOpts{:});

% Extract neural trajectories for original, unsegmented trials
% using learned parameters
[seqTrain, LLtrain] = exactInferenceWithLL_dlag(seqTrain, estParams, ...
                          'getLL', true);

% ==================================
% Assess generalization performance
% ==================================
if ~isempty(seqTest)
  %%% Conditional statistics (leave out select populations)
    % Pairwise regression on test data
    [~, MSE_reg, MSEorth_reg, R2_reg, R2orth_reg] ...
        = pairwise_regress_dlag(seqTest, estParams, 1:xDim_across, rGroups);
    
  %%% Joint statistics (evaluated on all populations jointly)
    % Log-likelihood and reconstruction error on test data
    % Including all types of latents (within and across)
    [~, LLtest.joint, ...
     MSE_denoise.joint, ~, ...
     R2_denoise.joint, ~] = denoise_dlag(seqTest, estParams);
 
  %%% Marginal statistics (evaluated on each population individually)
    % Partition input data and params according to observation groups.
    groupSeq = partitionObs(seqTest, yDims);
    groupParams = partitionParams_dlag(estParams);
    % Log-likelihood and reconstruction error on test data
    for groupIdx = 1:numGroups
        % Including all types of latents (within and across) 
        [~, LLtest.indiv(groupIdx), ...
         MSE_denoise.indiv(groupIdx), ...
         MSEorth_denoise.indiv{groupIdx}, ...
         R2_denoise.indiv(groupIdx), ...
         R2orth_denoise.indiv{groupIdx}] = denoise_dlag(groupSeq{groupIdx}, groupParams{groupIdx});
    end
end

% ===========
% Save results
% ===========

vars = who;
if ~parallelize
      fprintf('Saving %s...\n', fname);
end

% Remove redundant variables before saving.
vars = vars(~ismember(vars, {'varargin'}));
% Remove other unwanted variables
vars = vars(~ismember(vars, {'groupSeq', 'groupParams', 'numGroups', ...
                             'groupIdx'}));
% If saveData is true, then we'll write seqTrain, seqTest, seqTrainCut. 
% Otherwise, we won't.
if ~saveData
    vars = vars(~ismember(vars, {'seqTrain', 'seqTest', 'seqTrainCut'}));
end
% If previous model parameters were used for initialization, then startTau,
% etc. are not relevant
if isequal(init_method, 'params')
    vars = vars(~ismember(vars, {'startTau', 'startEps', 'startDelay'}));
end
save(fname,vars{:});
