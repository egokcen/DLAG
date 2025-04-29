function [res] = learnGPparams_plusDelays(seq, params, varargin)
%
% learnGPparams_plusDelays
%
% Description: Updates parameters of GP state model given inferred
%              latent states, jointly with time delays.
%              NOTE: Learning GP noise variance (eps) is currently
%                    unsupported.
%
% Arguments:
%
%     Required:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
% 
%     params -- Structure containing DLAG ACROSS-GROUP model parameters. 
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
%     18 Mar 2020 -- Initial full revision. 
%     19 Feb 2023 -- Added spectral Gaussian compatibility.

MAXITERS  = -10; % for minimize.m
verbose   = false;
assignopts(who, varargin);

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);

switch params.covType
    case 'rbf'
        % If there's more than one type of parameter, put them in the
        % second row of oldParams.
        oldParams = [params.gamma; params.DelayMatrix];
        fname     = 'grad_rbf_delay';
    case 'sg'
        % If there's more than one type of parameter, put them in the
        % second row of oldParams.
        oldParams = [params.gamma; params.DelayMatrix; params.nu];
        fname     = 'grad_sg_delay';
end

precomp = makePrecomp_timescales(seq,params);

% Loop once for each state dimension (each GP)
gamma = zeros(1,xDim);
DelayMatrix = zeros(numGroups, xDim);
nu = zeros(1,xDim);
for i = 1:xDim
    const = [];
    
    if ~params.notes.learnGPNoise
        const.eps = params.eps(i);
    end
    
    switch fname                
        case 'grad_rbf_delay'
            % Change of variables for constrained optimization
            init_gam = log(oldParams(1,i));
            % We don't include delays to the first group in the optimization
            init_delay = reshape(oldParams(3:end,i), numGroups-1, 1);
            init_delay = -log(2*params.maxDelay./(init_delay + params.maxDelay) - 1);
            init_p = [init_gam; init_delay];
            
        case 'grad_sg_delay'
            % Change of variables for constrained optimization
            init_gam = log(oldParams(1,i));
            % We don't include delays to the first group in the optimization
            init_delay = reshape(oldParams(3:(end-1),i), numGroups-1, 1);
            init_delay = -log(2*params.maxDelay./(init_delay + params.maxDelay) - 1);
            init_nu = oldParams(end,i);
            init_p = [init_gam; init_delay; init_nu];
    end   
    
    % This does the heavy lifting
    [res_p, fX, res_iters] =...
        minimize(init_p, fname, MAXITERS, precomp(i), const);
    
    switch params.covType
        case 'rbf'
            switch fname                
                case 'grad_rbf_delay'
                    gamma(i) = exp(res_p(1));
                    DelayMatrix(2:end,i) = 2*params.maxDelay./(1+exp(-res_p(2:end))) - params.maxDelay;
            end      
            
        case 'sg'               
            gamma(i) = exp(res_p(1));
            DelayMatrix(2:end,i) = 2*params.maxDelay./(1+exp(-res_p((1:numGroups-1)+1))) - params.maxDelay;
            nu(i) = res_p(end);
    end    
    
    if verbose
        fprintf('\nConverged p; xDim:%d, p:%s', i, mat2str(res_p, 3));
    end
end

res.DelayMatrix = DelayMatrix;
res.gamma = gamma;
if isequal(params.covType, 'sg')
    res.nu = nu; 
end
res.res_iters = res_iters; % Number of optimization iterations
res.fX = fX;               % Value of optimization objective function

end

