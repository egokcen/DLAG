function res = learnGPparams_freq(seq, Xspec, params, varargin)
%
% res = learnGPparams_freq(seq, Xspec, params, ...)
% 
% Description: Update parameters of GP state model using a frequency domain
%              approximation.
%
% Arguments:
%
%     Required:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId         -- unique trial identifier
%                     T (1 x 1)       -- number of timesteps
%                     yfft (yDim x T) -- unitary FFT neural data
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
%     verbose  -- logical that specifies whether to display status messages
%                 (default: false)
%
% Outputs:
%
%     res -- Structure containing the updated GP state model parameters.
%            Also includes the number of iterations and value of the cost 
%            function after updating these parameters via gradient descent.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     21 Jul 2023 -- Initial full revision.

  MAXITERS  = -10; % for minimize.m
  verbose   = false;
  assignopts(who, varargin);

  switch params.covType
    case 'rbf'
      % If there's more than one type of parameter, put them in the
      % second row of oldParams.
      oldParams = params.gamma;
      fname     = 'grad_rbf_freq';
    case 'sg'
      oldParams = [params.gamma; params.nu];
      fname     = 'grad_sg_freq';
  end
  if params.notes.learnGPNoise
    oldParams = [oldParams; params.eps];
    fname     = [fname '_noise'];
  end

  xDim    = size(oldParams, 2);
  precomp = makePrecomp_freq(seq, Xspec, xDim);
  
  % Loop once for each state dimension (each GP)
  for i = 1:xDim
    const = [];
    if ~params.notes.learnGPNoise  
      const.eps = params.eps(i);     
    end

    switch params.covType              
        case 'rbf'
            % Change of variables for constrained optimization
            initp = log(oldParams(1,i));
            
        case 'sg'
            % Change of variables for constrained optimization
            init_gam = log(oldParams(1,i));
            init_nu = oldParams(2,i);
            initp = [init_gam; init_nu];
    end

    % This does the heavy lifting
    [res_p, res_f, res_iters] =...
    minimize(initp, fname, MAXITERS, precomp(i), const);
    
    switch params.covType
      case 'rbf'
        res.gamma(i) = exp(res_p(1));
      case 'sg'
        res.gamma(i) = exp(res_p(1));
        res.nu(i) = res_p(2);
    end    
    if params.notes.learnGPNoise  
      res.eps(i) = exp(res_p(end));
    end
      
    if verbose
      fprintf('\nConverged p; xDim:%d, p:%s', i, mat2str(res_p, 3));
    end
  end
