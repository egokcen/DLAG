function bootParams = bootstrapGPparams(seq, params, binWidth, numBootstrap, varargin)
%
% bootParams = bootstrapGPparams(seq, params, binWidth, numBootstrap, ...)
%
% Description: Get an estimate of the uncertainty of DLAG Gaussian process 
%              timescales and delays using bootstrap samples.
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
%     binWidth         -- float; bin width or sample period, in units of
%                         time
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
%     segLength        -- int; length of segments to extract, in number of 
%                         timesteps. If infinite, entire trials are 
%                         extracted, i.e., no segmenting. (default: Inf)
%     tolLL            -- float; stopping criterion for EM (based on LL) 
%                         (default: 1e-4). 
%                         Note: The default here is larger than in 
%                         em_dlag.m, since it's assumed that an 
%                         already-fitted model is used for initialization.
%     maxIters         -- int; number of EM iterations to run 
%                         (default: 1e4)
%
% Outputs:
%
%     bootParams  -- Structure containing bootstrapped DLAG GP parameters.
%                    Contains the fields
%
%                    tau_across.upper -- (1 x xDim_across) array; 
%                                         confidence interval upper bound
%                                         on across-group GP timescales  
%                    tau_across.lower -- (1 x xDim_across) array; 
%                                         confidence interval lower bound
%                                         on across-group GP timescales 
%                    tau_within.upper -- (1 x numGroups) cell array; 
%                                        confidence interval upper bound
%                                        on GP timescales for each group
%                    tau_within.lower -- (1 x numGroups) cell array; 
%                                        confidence interval lower bound
%                                        on GP timescales for each group
%                    DelayMatrix.upper -- (numGroups x xDim_across) array;
%                                        confidence interval upper bound on 
%                                        the delay matrix, converted to 
%                                        units of time.
%                    DelayMatrix.lower -- (numGroups x xDim_across) array;
%                                        confidence interval lower bound on 
%                                        the delay matrix, converted to 
%                                        units of time.
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 May 2020 -- Initial full revision.
%     23 May 2020 -- Added parallel options.
%     26 May 2020 -- Added options for cut trials and EM convergence
%                    tolerance.
%     22 Sep 2020 -- Fixed issue with DelayMatrix shape for models with one
%                    across-group dimension.
%     20 Mar 2021 -- Added 'maxIters' option.

alpha       = 0.05;
parallelize = false;
numWorkers  = 4;
segLength   = Inf;
tolLL       = 1e-4;
maxIters    = 1e4;
extraOpts   = assignopts(who, varargin);

% Constants
learnObs = false;
verbose = false;
numGroups = length(params.yDims);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;

% Store results across all bootstrap samples
tau_across = nan(numBootstrap, xDim_across);
DelayMatrix = nan(numBootstrap, numGroups, xDim_across);
tau_within = cell(1,numGroups);
for groupIdx = 1:numGroups 
    tau_within{groupIdx} = nan(numBootstrap,xDim_within(groupIdx));
end

% For compute efficiency, train on equal-length segments of trials
seqCut = cutTrials(seq, 'segLength', segLength);
if isempty(seqCut)
    fprintf('WARNING: no segments extracted for re-fitting.  Defaulting to segLength=Inf.\n');
    seqCut = cutTrials(seq, 'segLength', Inf);
end
N = length(seqCut);

% Re-fit GP parameters on bootstrap samples, keeping the remaining
% parameters fixed.
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
    newParams = cell(1,numBootstrap);
    
    % Re-fit parameters
    parfor bootIdx = 1:numBootstrap
        fprintf('Re-fitting GP parameters on bootstrap sample %d of %d\n', ...
                bootIdx, numBootstrap);
        seqBoot = seqCut(bootSamples{bootIdx});
        [newParams{bootIdx}, ~, ~, ~, ~, ~, ~, ~, ~] = em_dlag(params, seqBoot, ...
                                                               'maxIters', maxIters, ...
                                                               'tolLL', tolLL, ...
                                                               'learnObs', learnObs, ...
                                                               'verbose', verbose);
        % Store results for the current bootstrap sample
        newParams{bootIdx} = getGPparams_dlag(newParams{bootIdx}, binWidth);
    end
    
    % Store results for the all bootstrap samples
    for bootIdx = 1:numBootstrap
        tau_across(bootIdx,:) = newParams{bootIdx}.tau_across;
        DelayMatrix(bootIdx,:,:) = newParams{bootIdx}.DelayMatrix;
        for groupIdx = 1:numGroups
            tau_within{groupIdx}(bootIdx,:) = newParams{bootIdx}.tau_within{groupIdx};
        end
    end
else
    for bootIdx = 1:numBootstrap
        fprintf('Re-fitting GP parameters on bootstrap sample %d of %d\n', ...
                bootIdx, numBootstrap);
        % Draw bootstrap sample
        bootSamples = datasample(1:N, N, 'Replace', true);
        seqBoot = seqCut(bootSamples);

        [newParams, ~, ~, ~, ~, ~, ~, ~, ~] = em_dlag(params, seqBoot, ...
                                                      'maxIters', maxIters, ...
                                                      'tolLL', tolLL, ...
                                                      'learnObs', learnObs, ...
                                                      'verbose', verbose);
        % Store results for the current bootstrap sample
        newParams = getGPparams_dlag(newParams, binWidth);
        tau_across(bootIdx,:) = newParams.tau_across;
        DelayMatrix(bootIdx,:,:) = newParams.DelayMatrix;
        for groupIdx = 1:numGroups
            tau_within{groupIdx}(bootIdx,:) = newParams.tau_within{groupIdx};
        end
    end
end

% Compute confidence intervals. Fill out output structures

% Across-group timescales
ci = prctile(tau_across, 100.*[alpha/2 1-alpha/2], 1);
bootParams.tau_across.lower = ci(1,:);
bootParams.tau_across.upper = ci(2,:);

% Within-group timescales
bootParams.tau_within.lower = cell(1,numGroups);
bootParams.tau_within.upper = cell(1,numGroups);
for groupIdx = 1:numGroups
    ci = prctile(tau_within{groupIdx}, 100.*[alpha/2 1-alpha/2], 1);
    bootParams.tau_within.lower{groupIdx} = ci(1,:); 
    bootParams.tau_within.upper{groupIdx} = ci(2,:);
end

% Delays
ci = prctile(DelayMatrix, 100.*[alpha/2 1-alpha/2], 1);
bootParams.DelayMatrix.lower = squeeze(ci(1,:,:));
bootParams.DelayMatrix.upper = squeeze(ci(2,:,:));
% 'squeeze' incorrectly shapes DelayMatrix into a row vector for 1D models
if xDim_across == 1 
    bootParams.DelayMatrix.upper = reshape(bootParams.DelayMatrix.upper,numGroups,[]);
    bootParams.DelayMatrix.lower = reshape(bootParams.DelayMatrix.lower,numGroups,[]);
end

