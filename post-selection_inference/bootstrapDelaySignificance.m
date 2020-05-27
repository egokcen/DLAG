function sig = bootstrapDelaySignificance(seq, params, numBootstrap, varargin)
%
% sig = bootstrapDelaySignificance(seq, params, numBootstrap, ...)
%
% Description: Evaluate the statistical signficance of each across-group
%              delay. "Significance" here is defined as the relative 
%              decrease in performance (relative to the unaltered model)
%              that results from setting a delay to zero. Get an estimate
%              of the variability of the prominence metric using bootstrap 
%              samples.
%
% Arguments:
%
%     Required:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
%
%     params  -- Structure containing DLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
%                    eps_across   -- (1 x xDim_across) GP noise variances
%                    gamma_within -- (1 x numGroups) cell array; 
%                                    GP timescales for each group
%                    eps_within   -- (1 x numGroups) cell array;
%                                    GP noise variances for each group
%                    d            -- (yDim x 1) array; observation mean
%                    C            -- (yDim x (numGroups*xDim)) array;
%                                    mapping between low- and high-d spaces
%                    R            -- (yDim x yDim) array; observation noise
%                                    covariance 
%                    DelayMatrix  -- (numGroups x xDim_across) array;
%                                    delays from across-group latents to 
%                                    observed variables. NOTE: Delays are
%                                    reported as (real-valued) number of
%                                    time-steps.
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%     numBootstrap     -- int; number of bootstrap samples
%
%     Optional:
%
%     alpha            -- float (in range [0 1]); significance level for 
%                         bootstrap confidence intervals
%     parallelize      -- logical; Set to true to use Matlab's parfor 
%                         construct to parallelize using multiple cores. 
%                         (default: false)
%     numWorkers       -- int; Number of cores to use, if using the 
%                         parallelize option. (default: 4)
%
% Outputs:
%
%     sig  -- structure containing the following fields:
%             raw    -- (1 x xDim_across) array; significance of each delay
%                       evaluated on raw data, measured by decrease in 
%                       log-likelihood relative to the unaltered model.
%             upper  -- (1 x xDim_across) array; upper bound of bootstrap 
%                       confidence interval
%             lower  -- (1 x xDim_across) array; lower bound of bootstrap 
%                       confidence interval
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 May 2020 -- Initial full revision.

alpha       = 0.05;
parallelize = false;
numWorkers  = 4;
extraOpts   = assignopts(who, varargin);

N = length(seq);
xDim_across = params.xDim_across;

% Evaluate significance on raw data
sig.raw = evalDelaySignificance(seq, params);

% Store results across all bootstrap samples
sig_all = nan(numBootstrap, xDim_across);

% Evaluate signficance on bootstrap samples
if parallelize
    % Set up parallelization, if desired
    StartParPool(numWorkers);  % Helper function to set up parfor construct
    % Draw all bootstrap samples up front, since parfor doesn't handle
    % randomization well.
    bootSamples = cell(1,numBootstrap);
    for bootIdx = 1:numBootstrap
        bootSamples{bootIdx} = datasample(1:N, N, 'Replace', true);
    end
    % Evaluate signficance
    parfor bootIdx = 1:numBootstrap
        seqBoot = seq(bootSamples{bootIdx});
        sig_all(bootIdx,:) = evalDelaySignificance(seqBoot, params);
    end
else
    for bootIdx = 1:numBootstrap
        % Draw bootstrap sample
        bootSamples = datasample(1:N, N, 'Replace', true);
        seqBoot = seq(bootSamples);
        sig_all(bootIdx,:) = evalDelaySignificance(seqBoot, params);
    end
end

% Compute confidence intervals. Fill out output structure.
ci = prctile(sig_all, 100.*[alpha/2 1-alpha/2], 1);
sig.lower = ci(1,:);
sig.upper = ci(2,:);
