function [seq, Corth] = predictiveProjection_dlag(seq, params, pDim, rGroups)
%
% seq = predictiveProjection_dlag(seq, params, pDim, rGroups)
%
% Description: Project activity of rGroups(1) onto the subspace most 
%              predictive of rGroups(2).
%              NOTE: Within-group (target) latent activity is excluded.
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
%     pDim    -- int; number of top orthonormal latent coordinates to use 
%                for projection
%     rGroups -- (1 x 2) int array; Indices of the groups (in yDims) we wish
%                to regress on. These groups should correspond to the blocks
%                of the matrix C
%
% OUTPUTS:
%
%     seq     -- test data structure with new field xorth_predXX, where XX is
%                pDim.
%     Corth   -- (yDims(rGroups(2)) x pDim) array; Top pDim predictive
%                dimensions
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     26 Mar 2020 -- Initial full revision.

  yDims = params.yDims;
  yDim = sum(yDims);
  xDim_across = params.xDim_across;
  xDim_within = params.xDim_within;
  numGroups = length(yDims);
  
  obs_block_idxs = get_block_idxs(yDims); % Index rows of C
  lat_block_idxs = get_block_idxs(xDim_across + xDim_within); % Index cols of C

  numRGroups = length(rGroups); % Really just 2, for pairwise regression

  fn = sprintf('xorth_pred%02d', pDim);
  for n = 1:length(seq)
    seq(n).(fn) = nan(pDim, seq(n).T);
  end
    
  % Indices of the current source group (predictors)
  currSourceGroup = rGroups(1);
  currSourceBlock = obs_block_idxs{currSourceGroup};
  sourceIdxs = currSourceBlock(1):currSourceBlock(2);
    
  % Indices for delayed latents corresponding to source group
  currSourceLatBlock = lat_block_idxs{currSourceGroup};
  sourceLatIdxs = currSourceLatBlock(1):currSourceLatBlock(2);
    
  % Indices of the current target group (responses)
  currTargetGroup = rGroups(2);
  currTargetBlock = obs_block_idxs{currTargetGroup};
  targetIdxs = currTargetBlock(1):currTargetBlock(2);
    
  % Indices for delayed latents corresponding to target group
  currTargetLatBlock = lat_block_idxs{currTargetGroup};
  targetLatIdxs = currTargetLatBlock(1):currTargetLatBlock(2);
  targetLatIdxs = targetLatIdxs(1:xDim_across); % Keep only across-group latents
  
  for n = 1:length(seq)
    seqReg(n).T = seq(n).T;
    seqReg(n).y = seq(n).y(sourceIdxs,:);
  end
  paramsReg   = params;
  paramsReg.C = params.C(sourceIdxs,:);
  paramsReg.d = params.d(sourceIdxs);
  paramsReg.R = params.R(sourceIdxs,sourceIdxs);
    
  % Infer latents given only the source group
  % Note that we're inferring within-group latents here, too, but we'll throw
  % away those signals in the projection step.
  seqReg = exactInferenceWithLL_dlag(seqReg, paramsReg, 'getLL', false);
    
  % Orthonormalize with respect to the *target* group
  Ctarget = params.C(targetIdxs,targetLatIdxs);
  Xtemp = [seqReg.xsm];
  [Xorth, Corth]   = orthogonalize(Xtemp(targetLatIdxs,:), Ctarget);
  seqReg           = segmentByTrial(seqReg, Xorth, 'xorth');  
  
  for n = 1:length(seq)
    seq(n).(fn) = seqReg(n).xorth(1:pDim,:);
  end
  
  Corth = Corth(:,1:pDim);
