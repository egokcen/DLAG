function [postparams, postsig, postIdxs] = postprocessDelays(params, sig, varargin)
%
% [postparams, postsig, postIdxs] = postprocessDelays(params, sig, ...)
%
% Description: Remove or set insignificant delays (determined by the 
%              bootstrapped results in sig) to 0. Sort latents according to
%              delay value, if desired.
%
% Arguments:
%
%     Required:
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
%                    xDim_within  -- (1 x numGroups) array; number of
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%     sig  -- structure containing the following fields:
%             raw    -- (1 x xDim_across) array; significance of each delay
%                       evaluated on raw data, measured by decrease in 
%                       log-likelihood relative to the unaltered model.
%             upper  -- (1 x xDim_across) array; upper bound of bootstrap 
%                       confidence interval
%             lower  -- (1 x xDim_across) array; lower bound of bootstrap 
%                       confidence interval
%
%     Optional:
%
%     remove  -- logical; set to true to remove across-group latents with 
%                delays that do not significantly differ from 0. Otherwise,
%                those delays will be set to 0, but kept among the
%                parameters (default: false)
%     sorted  -- logical; set to true to sort latents in descending order of
%                delay. rGroups might also need to be specified to indicate
%                delay directionality. (default: false)
%     rGroups -- (1 x 2) array; Indexes of groups to get relative
%                delays. rGroups(1) is the reference group (default: [1 2])
%
% Outputs:
%
%     postparams  -- Same format as params, modified according to delay 
%                    significance and sorting. 
%     postsig     -- Same format as sig, modified according to delay
%                    significance and sorting.
%     postIdxs    -- (1 x xDim_across) array; Indexes into original
%                    params/sig structure, prior to postprocessing. For
%                    example,
%                        postsig.lower = sig.lower(postIdxs);
%                    Useful if other post-selection inference structures,
%                    such as prominence results, need to be modified
%                    outside of this function.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     14 Sep 2020 -- Initial full revision.

remove = false;
sorted = false;
rGroups = [1 2];
extraOpts = assignopts(who, varargin);
numGroups = length(params.gamma_within);

% When postprocessing delays, we'll always leave within-group latents
% untouched.
withinIdxs = cell(1,numGroups);
for groupIdx = 1:numGroups
    withinIdxs{groupIdx} = 1:params.xDim_within(groupIdx); 
end

% Initialize output structures
postparams = params;
postsig = sig;
postIdxs = 1:params.xDim_across;

if sorted
    % Sort latents in descending delay order.
    % For sorting, take only the delays associated with the groups in rGroups
    delays = params.DelayMatrix(rGroups(2),:) - params.DelayMatrix(rGroups(1),:);
    [~, sortIdxs] = sort(delays, 'descend');
    postparams = getSubsetParams_dlag(postparams, sortIdxs, withinIdxs);
    postsig = getSubsetDelaySignificance(postsig, sortIdxs);
    postIdxs = sortIdxs;
end

if remove
    % Remove across-group latents associated with insignificant delays
    keptIdxs = find(postsig.lower > 0);
    postparams = getSubsetParams_dlag(postparams, keptIdxs, withinIdxs);
    postsig = getSubsetDelaySignificance(postsig, keptIdxs);
    postIdxs = postIdxs(keptIdxs);
else
    % Set insignificant delays of to zero
    postparams.DelayMatrix(:,(postsig.lower <= 0)) = 0;
end