function res = learnGPparams(seq, params, varargin)
% Updates parameters of GP state model given neural trajectories.
%
% INPUTS:
%
% seq         - data structure containing neural trajectories
% params      - current GP state model parameters, which gives starting point
%               for gradient optimization
%
% OUTPUT:
%
% res         - updated GP state model parameters
%
% OPTIONAL ARGUMENTS:
%
% MAXITERS    - maximum number of line searches (if >0), maximum number 
%               of function evaluations (if <0), for minimize.m (default:-8)
% verbose     - logical that specifies whether to display status messages
%               (default: false)
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
%     Modified from original code by Byron Yu and John Cunningham (2009).
% 
% Revision history:
%     19 Feb 2023 -- Changed kernel options to squared exponential (rbf)
%                    and spectral Gaussian (sg).

  MAXITERS  = -10; % for minimize.m
  verbose   = false;
  assignopts(who, varargin);

  switch params.covType
    case 'rbf'
      % If there's more than one type of parameter, put them in the
      % second row of oldParams.
      oldParams = params.gamma;
      fname     = 'grad_rbf';
    case 'sg'
      oldParams = [params.gamma; params.nu];
      fname     = 'grad_sg';
  end
  if params.notes.learnGPNoise
    oldParams = [oldParams; params.eps];
    fname     = [fname '_noise'];
  end

  xDim    = size(oldParams, 2);
  precomp = makePrecomp(seq, xDim);
  
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
