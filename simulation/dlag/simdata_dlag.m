function [seq, params] = simdata_dlag(N, T, binWidth, yDims, ...
    xDim_across, xDim_within, snr, covType, param_lims, varargin)
% [seq, params] = simdata_dlag(...)
%
% Description: Randomly generate DLAG model parameters, within specified
%              constraints. Then, generate synthetic data according to that
%              DLAG generative model.
%
% Arguments:
%
%     Required:
%
%     N         -- int; number of sequences
%     T         -- int; number of samples per sequence
%     binWidth  -- float; intended spike count bin width or sample period 
%                  (in units of time). Assume uniform sampling.
%     yDims     -- (1 x numGroups) array; List of dimensionalities of
%                  observed data, [Y1, Y2,...]
%     xDim_across -- int; Number of across-group latents
%     xDim_within -- (1 x numGroups) array; Number of within-group latents
%                    for each group. 
%     snr      -- (1 x numGroups) array; List of signal-to-noise ratios,
%                 defined as trace(CC') / trace(R)
%     covType  -- string; Specify type of covariance function to be used.
%                 Supported options:
%                     'rbf' -- Radial basis or squared exponential function
%                     'sg'  -- Spectral Gaussian or Gauss-cosine function
%     param_lims -- structure with fields that depend on the covType:
%         tau_lim   -- (1 x 2) array; lower- and upper-bounds of GP 
%                      timescales (for covTypes 'rbf' and 'sg')
%         eps_lim   -- (1 x 2) array; lower- and upper-bounds of GP noise 
%                      variances (for covTypes 'rbf' and 'sg')
%         delay_lim -- (1 x 2) array; lower- and upper-bounds of delays, in
%                      units of time (for covTypes 'rbf' and 'sg')
%         nu_lim    -- (1 x 2) array; lower- and upper-bounds of center
%                      frequencies, in units of 1/time (same units as
%                      delays and timescales) (for covType 'sg')
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
%                                     nu_across.*binWidth 
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
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     07 Apr 2020 -- Initial full revision.
%     17 Apr 2020 -- Added 0-within-group dimension functionality
%     28 Jun 2020 -- Added 0-across-group dimension functionality
%     14 Oct 2020 -- Added option to manually specify GP parameters
%     09 Jan 2021 -- Overhauled to utilize more modularized / elegant
%                    subroutines
%     18 Feb 2023 -- Added spectral Gaussian compatibility.

latentfield = 'xsm';
obsfield = 'y';
assignopts(who, varargin);

% Generate DLAG model parameters
params = generate_params_dlag(yDims, xDim_across, xDim_within, binWidth, ...
                              snr, covType, param_lims);

% Generate latent sequences
seq = generate_latents_dlag(params, T, N, 'latentfield', latentfield);

% Generate observed sequences
seq = generate_obs_dlag(seq, params, 'latentfield', latentfield, ...
                        'obsfield', obsfield);