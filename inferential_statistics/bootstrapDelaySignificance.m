function sig = bootstrapDelaySignificance(seq, params, numBootstrap, varargin)
%
% sig = bootstrapDelaySignificance(seq, params, numBootstrap, ...)
%
% Description: Evaluate the statistical signficance of each across-group
%              delay using a non-parametric bootstrap procedure. 
%              "Significance" here is defined as the relative decrease in 
%              performance (relative to the unaltered model) that results
%              from setting a delay to zero. Return the proportion of 
%              bootstrap samples on which the unaltered model gave a lower
%              likelihood than the zero-delay model.
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
%     parallelize      -- logical; Set to true to use Matlab's parfor 
%                         construct to parallelize using multiple cores. 
%                         (default: false)
%     numWorkers       -- int; Number of cores to use, if using the 
%                         parallelize option. (default: 4)
%
% Outputs:
%
%     sig  -- (1 x xDim_across) array; sig(i) gives the proportion of
%             bootstrap samples on which the model with delay i set to zero
%             performed at least as well as the unaltered model.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 May 2020 -- Initial full revision.
%     19 Jul 2021 -- Function now outputs the proportion of bootstrap 
%                    samples in which the zero-delay model performs at
%                    least as well as the original model, rather than the 
%                    alpha-percentile of the bootstrap distribution.

parallelize = false;
numWorkers  = 4;
extraOpts   = assignopts(who, varargin);

N = length(seq);
xDim_across = params.xDim_across;

% Store the bootstrap distributions for each delay
sig_dists = nan(numBootstrap, xDim_across);

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
        sig_dists(bootIdx,:) = evalDelaySignificance(seqBoot, params);
    end
else
    for bootIdx = 1:numBootstrap
        % Draw bootstrap sample
        bootSamples = datasample(1:N, N, 'Replace', true);
        seqBoot = seq(bootSamples);
        sig_dists(bootIdx,:) = evalDelaySignificance(seqBoot, params);
    end
end

% For each delay, report the proportion of bootstrap samples in which the
% zero-delay model performs at least as well as the original model.
sig = nan(1,xDim_across);
for dimIdx = 1:xDim_across
    sig(dimIdx) = mean(sig_dists(:,dimIdx) <= 0);
end
