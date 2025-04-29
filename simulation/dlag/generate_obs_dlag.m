function seq = generate_obs_dlag(seq, params, varargin)
% 
% seq = generate_obs_dlag(seq, params,...)
%
% Description: Generate observations given a set of latent sequences, all
%              from the DLAG generative model.
%
% Arguments:
%
%     Required:
%
%     seq     -- structure whose nth entry (corresponding to the nth 
%                sequence) has fields
%                    trialId   -- unique trial (sequence) identifier  
%                    T (1 x 1) -- number of timesteps
%                    (latentfield) (p x T) -- latent sequence
%
%     params  -- Structure containing DLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in units of time are given by 
%                                    'binWidth ./ sqrt(gamma)'                                                    
%                    eps_across   -- (1 x xDim_across) GP noise variances
%                    gamma_within -- (1 x numGroups) cell array; 
%                                    GP timescales for each group
%                    eps_within   -- (1 x numGroups) cell array;
%                                    GP noise variances for each group
%                    if covType == 'sg'
%                        nu_across -- (1 x xDim_across) array; center
%                                     frequencies for spectral Gaussians;
%                                     convert to 1/time via 
%                                     nu_across./binWidth 
%                        nu_within -- (1 x numGroups) cell array; 
%                                     center frequencies for each group
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
%     latentfield -- string; Name of latent data field in seq 
%                    (default: 'xsm')
%     obsfield    -- string; Name of observation data field in seq 
%                    (default: 'y')
%
% Outputs:
%     seq -- structure whose nth entry (corresponding to the nth 
%            sequence) has fields
%                trialId   -- unique trial (sequence) identifier  
%                T (1 x 1) -- number of timesteps
%                (latentfield) (p x T) -- latent sequence
%                (obsfield) (q x T) -- observation sequence
%
%
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     09 Jan 2021 -- Initial full revision.
%     08 Mar 2021 -- Updated documentation.

latentfield = 'xsm';
obsfield = 'y';
assignopts(who, varargin);

N = length(seq);
yDim = sum(params.yDims);

for n = 1:N
    T = seq(n).T;
    X = seq(n).(latentfield);
    ns = mvnrnd(zeros(1,yDim), params.R, T)';
    seq(n).(obsfield) = params.C * X + repmat(params.d,1,T) + ns;
end
