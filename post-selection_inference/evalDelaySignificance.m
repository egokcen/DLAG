function sig = evalDelaySignificance(seq, params)
%
% sig = evalDelaySignificance(seq, params)
%
% Description: Evaluate the statistical signficance of each across-group
%              delay. "Significance" here is defined as the relative 
%              decrease in performance (relative to the unaltered model)
%              that results from setting a delay to zero.
%
% Arguments:
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
%
% Outputs:
%
%     sig  -- (1 x xDim_across) array; sig(i) contains the significance 
%             of delays for across-group dimension i, measured by decrease 
%             in log-likelihood relative to the unaltered model.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 May 2020 -- Initial full revision.

xDim_across = params.xDim_across;

% Initialize output structure
sig = nan(1,xDim_across);

% Evaluate the likelihood of the full model, to provide a baseline.
[~, LL_full] = exactInferenceWithLL_dlag(seq, params, 'getLL', true);

% Evaluate delay signficance
for dimIdx = 1:xDim_across
    % Zero-out the current set of delays.
    currParams = params;
    currParams.DelayMatrix(:,dimIdx) = 0;
    % Evaluate performance
    [~, LL_zero] = exactInferenceWithLL_dlag(seq, currParams, 'getLL', true);
    % Store results
    sig(dimIdx) = LL_full - LL_zero;
end