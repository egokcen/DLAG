function [postparams, postprom, postIdxs_across, postIdxs_within] = postprocessProminence(params, prom, varargin)
%
% [postparams, postprom, postIdxs_across, postIdxs_within] = postprocessProminence(params, prom, ...)
%
% Description: Remove insignificant latent variables, according to overall
%              prominence. Sort latents according to prominence value, if 
%              desired.
%
%              NOTE: In general, prominence can be evaluated jointly on
%                    groups of latent variables, but this function assumes
%                    that prominence was evaluated for each latent variable
%                    individually.
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
%                 dimGroups -- (1 x numDimGroups) cell array; each element 
%                              contains an array of across-group latent 
%                              dimensions whose joint prominence was 
%                              evaluated.
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
%                 dimGroups -- (1 x numGroups) cell array; each element 
%                              contains a cell array of the same format as
%                              across.dimGroups, corresponding to groups of 
%                              within-group dimensions whose joint 
%                              prominence was evaluated.
%
%     Optional:
%
%     remove  -- logical; set to true to remove latents with prominence 
%                that does not significantly differ from 0.
%                (default: true)
%     sorted  -- logical; set to true to sort latents in descending order 
%                of prominence (default: true)
%     metric  -- string; Specify which measure of prominence to use: 'LL'
%                or 'VE' (default: 'LL')
%
% Outputs:
%
%     postparams  -- Same format as params, modified according to 
%                    overall prominence and sorting. 
%     postprom    -- Same format as prom, modified according to prominence
%                    sorting.
%     postIdxs_across -- (1 x xDim_across) array; Indexes into original
%                    params/prom structure, prior to postprocessing. For
%                    example,
%                        postprom.across.LL.lower = prom.across.LL.lower(postIdxs_across);
%                    Useful if other post-selection inference structures,
%                    such as group prominence results, need to be modified
%                    outside of this function.
%     postIdxs_within -- (1 x numGroups) cell array; postIdxs_within{i} 
%                    indexes latents in kept in group i into the orginal
%                    params/prom structure.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     14 Sep 2020 -- Initial full revision.

remove = true;
sorted = true;
metric = 'LL';
extraOpts = assignopts(who, varargin);
numGroups = length(params.gamma_within);

% Initialize output structures
postparams = params;
postprom = prom;
postIdxs_across = 1:params.xDim_across;
postIdxs_within = cell(1,numGroups);
for groupIdx = 1:numGroups
    postIdxs_within{groupIdx} = 1:params.xDim_within(groupIdx); 
end

if sorted
    % Sort latents in descending order of prominence.
    [~, sortIdxs_across] = sort(postprom.across.(metric).raw, 'descend');
    sortIdxs_within = cell(1,numGroups);
    for groupIdx = 1:numGroups
        [~, sortIdxs_within{groupIdx}] = sort(postprom.within.(metric).raw{groupIdx}, 'descend');
    end
    postparams = getSubsetParams_dlag(postparams, sortIdxs_across, sortIdxs_within);
    postprom = getSubsetProminence(postprom, sortIdxs_across, sortIdxs_within);
    postIdxs_across = sortIdxs_across;
    postIdxs_within = sortIdxs_within;
end

if remove
    % Remove latents with insignificant overall prominence
    keptIdxs_across = find(postprom.across.(metric).lower > 0);
    keptIdxs_within = cell(1,numGroups);
    for groupIdx = 1:numGroups
        keptIdxs_within{groupIdx} = find(postprom.within.(metric).lower{groupIdx} > 0);
    end
    postparams = getSubsetParams_dlag(postparams, keptIdxs_across, keptIdxs_within);
    postprom = getSubsetProminence(postprom, keptIdxs_across, keptIdxs_within);
    postIdxs_across = postIdxs_across(keptIdxs_across);
    for groupIdx = 1:numGroups
        postIdxs_within{groupIdx} = postIdxs_within{groupIdx}(keptIdxs_within{groupIdx}); 
    end
end