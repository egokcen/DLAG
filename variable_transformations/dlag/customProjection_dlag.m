function seq = customProjection_dlag(seq, params, U)
%
% seq = customProjection_dlag(seq, params, U)
%
% Description: Project latent trajectories inferred by a DLAG model onto
%              the modes, specified in U.
%
% Arguments:
%
%     Required:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId        -- unique trial identifier
%                     T (1 x 1)      -- number of timesteps
%                     xsm (xDim x T) -- latent trajectories
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
%
%     U       -- (1 x 2) cell array; U{i} -- (yDims(i) x r(i)) array; a 
%                basis for a subspace in group i, of dimension r(i).
%
% Outputs:
%
%     seq     -- input data structure with new field
%                  xproj -- (r x T) array; trajectories (of all groups)
%                           projected onto specified modes.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     22 Feb 2023 -- Initial full revision.

numGroups = length(params.yDims);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
xDim_total = xDim_across + xDim_within;
N = length(seq);
  
% Initialize output structure
fn = 'xproj';
for n = 1:N
    seq(n).(fn) = [];
end

% Get appropriate latent trajectories
groupSeq = partitionObs(seq, xDim_total, 'datafield', 'xsm');
% Project latents for each group separately
groupParams = partitionParams_dlag(params);
for groupIdx = 1:numGroups
    
    C = groupParams{groupIdx}.C;
    X = [groupSeq{groupIdx}.xsm];
    
    % Project latents onto modes
    Xproj = U{groupIdx}' * C *X;
    
    % Add projected latents to output structure
    seqTemp = segmentByTrial(seq,Xproj,fn);
    for n = 1:N
        seq(n).(fn) = [seq(n).(fn); seqTemp(n).(fn)];
    end

end