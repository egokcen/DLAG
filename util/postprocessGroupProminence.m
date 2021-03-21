function [postparams, postprom, postIdxs_across, postIdxs_within] = postprocessGroupProminence(params, prom, varargin)
%
% [postparams, postprom, postIdxs_across, postIdxs_within] = postprocessGroupProminence(params, prom, ...)
%
% Description: Remove insignificant latent variables, according to
%              individual group prominence. Sort latents according to 
%              prominence value, if desired.
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
%     prom -- (1 x numGroups) cell array; inprom{i} contains a structure
%             with the bootstrapped prominence of each group of latent 
%             variables for observation group i. See bootstrapProminence 
%             for details on the format of these structures.
%
%     Optional:
%
%     remove    -- logical; set to true to remove latents with prominence 
%                  that does not significantly differ from 0.
%                  (default: true)
%     sorted    -- logical; set to true to sort latents in descending order 
%                  of prominence (default: true)
%     sortGroup -- int; if 'sorted' is true, specifies the group with
%                  respect to which across-group latents are being sorted.
%     metric    -- string; Specify which measure of prominence to use: 'LL'
%                  or 'VE' (default: 'LL')
%
% Outputs:
%
%     postparams  -- Same format as params, modified according to 
%                    individual group prominence and sorting. 
%     postprom    -- Same format as prom, modified according to prominence
%                    sorting.
%     postIdxs_across -- (1 x xDim_across) array; Indexes into original
%                    params/prom structure, prior to postprocessing. For
%                    example,
%                        postprom{i}.across.LL.lower = prom{i}.across.LL.lower(postIdxs_across);
%                    Useful if other post-selection inference structures,
%                    such as delay significance results, need to be modified
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
sortGroup = 1;
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
    % Sort latents in descending order of prominence in sortGroup.
    [~, sortIdxs_across] = sort(postprom{sortGroup}.across.(metric).raw, 'descend');
    sortIdxs_within = cell(1,numGroups);
    for groupIdx = 1:numGroups
        [~, sortIdxs_within{groupIdx}] = sort(postprom{groupIdx}.within.(metric).raw{1}, 'descend');
    end
    postparams = getSubsetParams_dlag(postparams, sortIdxs_across, sortIdxs_within);
    postprom = getSubsetGroupProminence(postprom, sortIdxs_across, sortIdxs_within);
    postIdxs_across = sortIdxs_across;
    postIdxs_within = sortIdxs_within;
end

if remove
    % Remove latents with insignificant prominence in at least one group
    keptAcross = ones(1,postparams.xDim_across);
    keptIdxs_within = cell(1,numGroups);
    for groupIdx = 1:numGroups
        keptAcross(postprom{groupIdx}.across.(metric).lower <= 0) = 0;
        keptIdxs_within{groupIdx} = find(postprom{groupIdx}.within.(metric).lower{1} > 0);
    end
    keptIdxs_across = find(keptAcross > 0);
    postparams = getSubsetParams_dlag(postparams, keptIdxs_across, keptIdxs_within);
    postprom = getSubsetGroupProminence(postprom, keptIdxs_across, keptIdxs_within);
    postIdxs_across = postIdxs_across(keptIdxs_across);
    for groupIdx = 1:numGroups
        postIdxs_within{groupIdx} = postIdxs_within{groupIdx}(keptIdxs_within{groupIdx}); 
    end
end