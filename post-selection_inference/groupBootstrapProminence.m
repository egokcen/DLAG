function prom = groupBootstrapProminence(seq, params, dimGroups_across, dimGroups_within, numBootstrap, varargin)
%
% prom = groupBootstrapProminence(seq, params, dimGroups_across, dimGroups_within, numBootstrap, ...)
%
% Description: Within each group of observations, evaluate the "prominence"
%              of each group of latents given in dimGroups_across and 
%              dimGroups_within. "Prominence" here is defined as the 
%              relative decrease in performance that results from removing 
%              a group of latents from the full model. Get an estimate of 
%              the variability of the prominence metric using bootstrap 
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
%     dimGroups_across -- (1 x numDimGroups) cell array; each element 
%                         contains an array of across-group latent dimensions
%                         whose joint prominence will be evaluated.
%     dimGroups_within -- (1 x numGroups) cell array; each element 
%                         contains a cell array of the same format as
%                         dimGroups_across, corresponding to groups of 
%                         within-group dimensions whose joint prominence
%                         will be evaluated.
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
%     prom -- (1 x numGroups) cell array; prom{i} contains a structure with
%             the bootstrapped prominence of each group of latent variables
%             for observation group i. See bootstrapProminence for details
%             on the format of these structures.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 May 2020 -- Initial full revision.
%     23 May 2020 -- Added parallel options.

alpha       = 0.05;
parallelize = false;
numWorkers  = 4;
extraOpts   = assignopts(who, varargin);

numGroups = length(params.yDims);

% Partition input data and params according to observation groups.
groupSeq = partitionObs(seq, params.yDims);
groupParams = partitionParams_dlag(params);
prom = cell(1,numGroups);
for groupIdx = 1:numGroups
    prom{groupIdx} = bootstrapProminence(groupSeq{groupIdx}, ...
                                         groupParams{groupIdx}, ...
                                         dimGroups_across, ...
                                         dimGroups_within(groupIdx), ...
                                         numBootstrap, ...
                                         'alpha', alpha, ...
                                         'parallelize', parallelize, ...
                                         'numWorkers', numWorkers);
end