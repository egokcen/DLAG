function prom = bootstrapProminence(seq, params, dimGroups_across, dimGroups_within, numBootstrap, varargin)
%
% prom = bootstrapProminence(seq, params, dimGroups_across, dimGroups_within, numBootstrap, ...)
%
% Description: Evaluate the "prominence" of each group of latents given in 
%              dimGroups_across and dimGroups_within. "Prominence" here
%              is defined as the relative decrease in performance that
%              results from removing a group of latents from the full
%              model. Get an estimate of the variability of the prominence
%              metric using bootstrap samples.
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
%     prom -- structure containing the following fields:
%             across -- structure with across-group prominence. Has fields
%                 LL.raw    -- (1 x numDimGroups) array; prominence of
%                              dimGroups_across(i) evaluated on raw data, 
%                              measured by decrease in log-likelihood 
%                              relative to the full model.
%                 LL.upper  -- (1 x numDimGroups) array; upper bound of
%                              bootstrap confidence interval
%                 LL.lower  -- (1 x numDimGroups) array; lower bound of
%                              bootstrap confidence interval
%                 VE.raw    -- (1 x numDimGroups) array;prominence of
%                              dimGroups_across(i) evaluated on raw data, 
%                              measured by normalized decrease in variance 
%                              explained relative to the full model.
%                 VE.upper  -- (1 x numDimGroups) array; upper bound of
%                              bootstrap confidence interval
%                 VE.lower  -- (1 x numDimGroups) array; lower bound of
%                              bootstrap confidence interval
%                 dimGroups -- (1x numDimGroups) cell array; Exactly
%                              dimGroups_across. Entries correspond to
%                              above fields.
%             within -- structure with within-group prominence. Has fields
%                 LL.raw    -- (1 x numGroups) cell array; prominence of 
%                              dimGroups_within{i} evaluated on raw data, 
%                              measured by decrease in log-likelihood 
%                              relative to the full model.
%                 LL.upper  -- (1 x numGroups) cell array; upper bound of
%                              bootstrap confidence interval
%                 LL.lower  -- (1 x numGroups) cell array; lower bound of
%                              bootstrap confidence interval
%                 VE.raw    -- (1 x numGroups) cell array; prominence of
%                              dimGroups_within{i} evaluated on raw data, 
%                              measured by normalized decrease in variance
%                              explained relative to the full model.
%                 VE.upper  -- (1 x numGroups) cell array; upper bound of
%                              bootstrap confidence interval
%                 VE.lower  -- (1 x numGroups) cell array; lower bound of
%                              bootstrap confidence interval
%                 dimGroups -- (1x numGroups) cell array; Exactly
%                              dimGroups_within. Entries correspond to
%                              above fields.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     15 May 2020 -- Initial full revision.
%     23 May 2020 -- Added parallel options.

alpha       = 0.05;
parallelize = false;
numWorkers  = 4;
extraOpts   = assignopts(who, varargin);

N = length(seq);
numGroups = length(params.yDims);

% Evaluate prominence on raw data
currProm = evalProminence(seq, params, dimGroups_across, dimGroups_within);
% Fill out output structures
prom.across.LL.raw = currProm.across.LL;
prom.across.VE.raw = currProm.across.VE;
prom.within.LL.raw = cell(1,numGroups);
prom.within.VE.raw = cell(1,numGroups);
for groupIdx = 1:numGroups
    prom.within.LL.raw{groupIdx} = currProm.within.LL{groupIdx};
    prom.within.VE.raw{groupIdx} = currProm.within.VE{groupIdx};
end

% Store results across all bootstrap samples
LL_across = nan(numBootstrap, length(dimGroups_across));
VE_across = nan(numBootstrap, length(dimGroups_across));
LL_within = cell(1,numGroups);
VE_within = cell(1,numGroups);
for groupIdx = 1:numGroups 
    LL_within{groupIdx} = nan(numBootstrap,length(dimGroups_within{groupIdx}));
    VE_within{groupIdx} = nan(numBootstrap,length(dimGroups_within{groupIdx}));
end

% Evaluate prominence on bootstrap samples
if parallelize
    % Set up parallelization, if desired
    StartParPool(numWorkers);  % Helper function to set up parfor construct
    % Draw all bootstrap samples up front, since parfor doesn't handle
    % randomization well.
    bootSamples = cell(1,numBootstrap);
    for bootIdx = 1:numBootstrap
        bootSamples{bootIdx} = datasample(1:N, N, 'Replace', true);
    end
    % Initialize parfor compatible output structure
    currProm = cell(1,numBootstrap);
    
    % Evaluate prominence
    parfor bootIdx = 1:numBootstrap
        seqBoot = seq(bootSamples{bootIdx});
        currProm{bootIdx} = evalProminence(seqBoot, params, dimGroups_across, dimGroups_within);
    end
    
    % Store results for the all bootstrap samples
    for bootIdx = 1:numBootstrap
        LL_across(bootIdx,:) = currProm{bootIdx}.across.LL;
        VE_across(bootIdx,:) = currProm{bootIdx}.across.VE;
        for groupIdx = 1:numGroups
            LL_within{groupIdx}(bootIdx,:) = currProm{bootIdx}.within.LL{groupIdx};
            VE_within{groupIdx}(bootIdx,:) = currProm{bootIdx}.within.VE{groupIdx};
        end 
    end
else
    for bootIdx = 1:numBootstrap
        % Draw bootstrap sample
        bootSamples = datasample(1:N, N, 'Replace', true);
        seqBoot = seq(bootSamples);
        currProm = evalProminence(seqBoot, params, dimGroups_across, dimGroups_within);

        % Store results for the current bootstrap sample
        LL_across(bootIdx,:) = currProm.across.LL;
        VE_across(bootIdx,:) = currProm.across.VE;
        for groupIdx = 1:numGroups
            LL_within{groupIdx}(bootIdx,:) = currProm.within.LL{groupIdx};
            VE_within{groupIdx}(bootIdx,:) = currProm.within.VE{groupIdx};
        end
    end
end

% Compute confidence intervals. Fill out output structure.
prom.across.dimGroups = dimGroups_across;
prom.within.dimGroups = dimGroups_within;

% LL
ci = prctile(LL_across, 100.*[alpha/2 1-alpha/2], 1);
prom.across.LL.lower = ci(1,:);
prom.across.LL.upper = ci(2,:);

prom.within.LL.lower = cell(1,numGroups);
prom.within.LL.upper = cell(1,numGroups);
for groupIdx = 1:numGroups
    ci = prctile(LL_within{groupIdx}, 100.*[alpha/2 1-alpha/2], 1);
    prom.within.LL.lower{groupIdx} = ci(1,:); 
    prom.within.LL.upper{groupIdx} = ci(2,:);
end

% VE
ci = prctile(VE_across, 100.*[alpha/2 1-alpha/2], 1);
prom.across.VE.lower = ci(1,:);
prom.across.VE.upper = ci(2,:);

prom.within.VE.lower = cell(1,numGroups);
prom.within.VE.upper = cell(1,numGroups);
for groupIdx = 1:numGroups
    ci = prctile(VE_within{groupIdx}, 100.*[alpha/2 1-alpha/2], 1);
    prom.within.VE.lower{groupIdx} = ci(1,:); 
    prom.within.VE.upper{groupIdx} = ci(2,:);
end
