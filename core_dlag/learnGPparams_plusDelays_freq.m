function [res] = learnGPparams_plusDelays_freq(seq, Xspec, params, varargin)
%
% learnGPparams_plusDelays_freq
%
% Description: Updates parameters of GP state model given inferred
%              latent states, jointly with time delays, using a
%              frequency domain approximation.
%
% Arguments:
%
%     Required:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId         -- unique trial identifier
%                     T (1 x 1)       -- number of timesteps
%                     yfft (yDim x T) -- unitary FFT of the neural data
%                     xfft (xDim x T) -- unitary FFT of the latent 
%                                        posterior mean at each frequency
%     Xspec     -- data structure whose jth entry, corresponding
%                  to a group of trials of the same length, has fields
%                      T       -- int; number of time steps for this
%                                 trial group
%                      Sx_post -- (xDim x xDim x T) array; posterior 
%                                 spectrum at each frequency
%                      NOTE: For DLAG, posterior covariance/spectra of X 
%                            are the same for trials of the same length.
% 
%     params -- Structure containing DLAG model parameters. 
%               Contains the fields
% 
%               covType -- string; type of GP covariance (e.g., 'rbf')
%               gamma   -- (1 x xDim) array; GP timescales
%                          in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
%               eps     -- (1 x xDim) GP noise variances
%               if covType == 'sg'
%                   nu -- (1 x xDim) array; center frequencies for spectral
%                         Gaussians; convert to 1/time via 
%                         nu_across./binWidth
%               d            -- (yDim x 1) array; observation mean
%               C            -- (yDim x (numGroups*xDim)) array;
%                               mapping between low- and high-d spaces
%               R            -- (yDim x yDim) array; observation noise
%                               covariance 
%               DelayMatrix  -- (numGroups x xDim_across) array;
%                               delays from across-group latents to 
%                               observed variables. NOTE: Delays are
%                               reported as (real-valued) number of
%                               time-steps.
%               xDim    -- int; number of across-group latent variables
%               yDims        -- (1 x numGroups) array; 
%                               dimensionalities of each observed group
%
%     Optional:
%
%     MAXITERS -- int; maximum number of line searches (if >0), 
%                 maximum number of function evaluations (if <0), 
%                 for minimize.m (default:-10)
%     learnDelays -- logical; set true to learn delay parameters;
%                    otherwise, delays will remain fixed at their initial
%                    value (default: true)
%
% Outputs:
%
%     res -- Structure containing the updated GP state model parameters.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     24 Jul 2023 -- Initial full revision.

MAXITERS  = -10; % for minimize.m
learnDelays = true;
assignopts(who, varargin);

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);

precomp = makePrecomp_freq(seq,Xspec,xDim);

switch params.covType
    case 'rbf'

        % ======================
        % Timescale parameters
        % ======================

        % Loop once for each state dimension (each GP)
        gamma = zeros(1,xDim);

        for i = 1:xDim
            const = [];

            if ~params.notes.learnGPNoise
                const.eps = params.eps(i);
            end

            init_gam = log(params.gamma(i));

            % This does the heavy lifting
            [res_gam,~,~] = minimize(init_gam, ...
                                     'grad_rbf_freq', ...
                                     MAXITERS, ...
                                     precomp(i), ...
                                     const);
            gamma(i) = exp(res_gam);

        end

        res.gamma = gamma;

    case 'sg'

        % ===========================================
        % Timescale and center frequency parameters
        % ===========================================

        % Loop once for each state dimension (each GP)
        gamma = zeros(1,xDim);
        nu = zeros(1,xDim);

        for i = 1:xDim
            const = [];

            if ~params.notes.learnGPNoise
                const.eps = params.eps(i);
            end

            init_gam = log(params.gamma(i));
            init_nu = params.nu(i);
            init_p = [init_gam; init_nu];

            % This does the heavy lifting
            [res_p,~,~] = minimize(init_p, ...
                                   'grad_sg_freq', ...
                                   MAXITERS, ...
                                   precomp(i), ...
                                   const);
            gamma(i) = exp(res_p(1));
            nu(i) = res_p(2);

        end

        res.gamma = gamma;
        res.nu = nu;
end

% ======================
% Delay parameters
% ======================

if learnDelays && numGroups > 1

    precomp = makePrecomp_delays_freq(seq,Xspec,params);

    DelayMatrix = zeros(numGroups, xDim);

    % We don't include delays to the first group in the optimization
    const = [];
    const.maxDelay = params.maxDelay;
    for groupIdx = 2:numGroups        

        init_delay = params.DelayMatrix(groupIdx,:).';
        init_delay = log(const.maxDelay + init_delay) ...
                   - log(const.maxDelay - init_delay);
        
        % This does the heavy lifting
        % This does the heavy lifting
        [res_delay,~,~] = minimize(init_delay, ...
                                    'grad_delay_freq', ...
                                    MAXITERS, ...
                                    precomp(groupIdx), ...
                                    const);
        
        DelayMatrix(groupIdx,:) = const.maxDelay.*tanh(res_delay./2);
    
    end

    res.DelayMatrix = DelayMatrix;
else
    res.DelayMatrix = params.DelayMatrix;
end
