function [seq, sortParams] = scaleByVarExp(seq, params, varexp, varargin)
%
% [seq, sortParams] = scaleByVarExp(seq, params, varexp, varargin)
%
% Description: Scale (and re-order) latent trajectories by their fraction 
%              of variance explained, for visualization purposes. 
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
%     varexp   -- (1 x numGroups) cell array; varexp{i} is an
%                 array with the variance explained by each
%                 individual latent variable.
%
%     Optional:
%
%     sortDims     -- logical; set to true to sort dimensions by variance
%                     explained. (default: true)
%     sortGroup    -- int; if sorting, choose which group to use as an
%                     anchor when sorting across-group dimensions.
%                     (default: 1)
%     numAcross    -- int; Keep at most the top numAcross across-group
%                     dimensions (default: xDim_across)
%     numWithin    -- int; Keep at most the top numWithin within-group 
%                     dimensions for each group. 
%                     (default: max(xDim_within))
%     indatafield  -- string; fieldname of input data (default: 'xsm')
%     outdatafield -- string; fieldname of output data (default: 'xve')
%
% Outputs:
%
%     seq     -- input data structure with new field
%                    (outdatafield) -- (xDim x T) array; scaled (and
%                                      re-ordered) latent trajectories
%     sortParams -- same structure as params, with dimensions re-sorted.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     19 Mar 2021 -- Initial full revision.

sortDims = true;
sortGroup = 1;
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
numAcross = xDim_across;
numWithin = max(xDim_within);
indatafield = 'xsm';
outdatafield = 'xve';
assignopts(who, varargin);

numGroups = length(params.yDims);
N = length(seq);

% Initialize output structures
sortParams = params;

% Scale latents by variance explained
varexp_all = [varexp{:}];
for n = 1:N
    seq(n).(outdatafield) = varexp_all' .* seq(n).(indatafield);
end

if sortDims
    % Initialize sort indices
    sortIdxs_across = 1:params.xDim_across;
    sortIdxs_within = cell(1,numGroups);
    for groupIdx = 1:numGroups
        sortIdxs_within{groupIdx} = 1:params.xDim_within(groupIdx); 
    end
    
    % Sort dimensions by variance explained. Sort within- and across-group
    % dimensions separately.
    
    % Sort all across-group dimensions according to the anchor group 
    [~, sortIdxs_across] = sort(varexp{sortGroup}(1:xDim_across),'descend');
    % Keep at most the top numAcross dimensions
    keptAcross = min([numAcross xDim_across]);
    sortIdxs_across = sortIdxs_across(1:keptAcross);
    sortIdxs_all = cell(1,numGroups);
    for groupIdx = 1:numGroups
        % Sort within-group dimensions for each group separately
        [~,sortIdxs_within{groupIdx}] = sort(varexp{groupIdx}(xDim_across+1:end),'descend');
        % Keep at most the top numWithin dimensions
        keptWithin = min([numWithin xDim_within(groupIdx)]);
        sortIdxs_within{groupIdx} = sortIdxs_within{groupIdx}(1:keptWithin);
        % Collect all sort indexes
        sortIdxs_all{groupIdx} = [sortIdxs_across xDim_across+sortIdxs_within{groupIdx}];
    end
    
    % Sort parameters
    sortParams = getSubsetParams_dlag(params, sortIdxs_across, sortIdxs_within);
    
    % Sort time courses
    groupSeq = partitionObs(seq,xDim_across+xDim_within,'datafield',outdatafield);
    for n = 1:N
        currSeq = [];
        for groupIdx = 1:numGroups
            currSeq = [currSeq; groupSeq{groupIdx}(n).(outdatafield)(sortIdxs_all{groupIdx},:)];
        end
        seq(n).(outdatafield) = currSeq;
    end
end