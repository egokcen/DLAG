function Ka = getZeroLagK_across(params, varargin)
%
% Ka = getZeroLagK_across(params, ...)
%
% Description: Given a pair of groups, extract the 0-lag portion of the 
%              across-group GP kernel matrix, Ka(t,t).
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
%
%     Optional:
%
%     groupIdxs -- (1 x 2) int array; Specify which pair of groups to
%                  analyze. Order doesn't matter. (default: [1 2])
%
% Outputs:
%
%     Ka -- (xDim_across x xDim_across) array; diagonal matrix with 0-lag
%           portion of the across-grop GP kernel matrix, i.e., the
%           cross-correlation between a pair of groups' latents when t1 =
%           t2 = t.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 Mar 2021 -- Initial full revision.

groupIdxs = [1 2];
assignopts(who,varargin);

% Get the 0-lag portion of the across-group GP kernel matrix, Ka(t,t)
xDim_across = params.xDim_across;
params_pair.xDim = xDim_across;
params_pair.yDims = params.yDims(groupIdxs);
params_pair.DelayMatrix = params.DelayMatrix(groupIdxs,:);
params_pair.gamma = params.gamma_across;
params_pair.eps = params.eps_across;
params_pair.covType = params.covType;
K_big = make_K_big_plusDelays(params_pair,1);
Ka = K_big(1:xDim_across, xDim_across+1:end);