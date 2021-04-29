function [seq, sortParams] = scaleByCorr(seq, params, lcorr, groupIdxs, varargin)
%
% [seq, sortParams] = scaleByCorr(seq, params, lcorr, groupIdxs, varargin)
%
% Description: Scale (and re-order) across-group trajectories by their 
%              of correlation, for visualization purposes. 
%
% Arguments:
%
%     Required:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId        -- unique trial identifier
%                     T (1 x 1)      -- number of timesteps
%                     (indatafield) (xDim x T) -- latent trajectories
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
%     lcorr   -- (1 x xDim_across) array; correlation within each 
%                across-group latent dimension
%     groupIdxs -- (1 x 2) array; Specify which groups were used to compute
%                  lcorr.
%
%     Optional:
%
%     sortDims     -- logical; set to true to sort dimensions by 
%                     correlation. (default: true)
%     numAcross    -- int; Keep at most the top numAcross across-group
%                     dimensions (default: xDim_across)
%     indatafield  -- string; fieldname of input data (default: 'xsm')
%     outdatafield -- string; fieldname of output data (default: 'xce')
%
% Outputs:
%
%     seq     -- input data structure with new field
%                    (outdatafield) -- (2*xDim_across x T) array; scaled 
%                                      (and re-ordered) latent trajectories
%     sortParams -- same structure as params, with dimensions re-sorted.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     19 Mar 2021 -- Initial full revision.
%     28 Apr 2021 -- Added exception handling for xDim_across = 0 case.

sortDims = true;
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
numAcross = xDim_across;
indatafield = 'xsm';
outdatafield = 'xce';
assignopts(who, varargin);

N = length(seq);
numGroups = length(params.yDims);

% Initialize output structure
sortParams = params;
for n = 1:N
    seq(n).(outdatafield) = []; 
end

if xDim_across <= 0
   fprintf('scaleByCorr: xDim_across = 0. Returning empty field seq.%s and sortParams unmodified.\n', outdatafield);
   return;
end

% Extract across-group latents, and keep only the pair of groups of interest
groupSeq = partitionObs(seq,xDim_across+xDim_within,'datafield',indatafield);
for groupIdx = groupIdxs
    [x_across, ~] = partitionLatents_meanOnly(groupSeq{groupIdx},xDim_across,xDim_within(groupIdx));
    for n = 1:N
        seq(n).(outdatafield) = [seq(n).(outdatafield); x_across(n).xsm];
    end
end
% Update other variables for smaller number of groups
numGroups = length(groupIdxs);
xDim_within = xDim_within(groupIdxs);

% Scale latents by correlation
for n = 1:N
    seq(n).(outdatafield) = repmat(lcorr,1,numGroups)' .* seq(n).(outdatafield);
end

if sortDims
    % Initialize sort indices
    sortIdxs_across = 1:params.xDim_across;
    
    % Sort dimensions by correlation 
    [~, sortIdxs_across] = sort(lcorr(1:xDim_across),'descend');
    % Keep at most the top numAcross dimensions
    keptAcross = min([numAcross xDim_across]);
    sortIdxs_across = sortIdxs_across(1:keptAcross);
    
    % Sort parameters
    sortParams = getSubsetParams_dlag(params,sortIdxs_across,cell(1,numGroups));
    
    % Sort time courses
    groupSeq = partitionObs(seq,repmat(xDim_across,1,numGroups),'datafield',outdatafield);
    for n = 1:N
        currSeq = [];
        for groupIdx = 1:numGroups
            currSeq = [currSeq; groupSeq{groupIdx}(n).(outdatafield)(sortIdxs_across,:)];
        end
        seq(n).(outdatafield) = currSeq;
    end
end