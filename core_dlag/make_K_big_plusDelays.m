function [K_big] = make_K_big_plusDelays(params,T)
%
% [K_big] = make_K_big_plusDelays(params, T)
%
% Description: Constructs full GP covariance matrix across all across-group
%              state dimensions and timesteps.
%
% Arguments:
%
%     params  -- Structure containing DLAG model parameters at which EM
%                algorithm is initialized. Contains (at minimum) the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma        -- (1 x xDim_across) array; GP timescales
%                                    in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
%                    eps          -- (1 x xDim_across) GP noise variances
%                    DelayMatrix  -- (numGroups x xDim_across) array;
%                                    delays from across-group latents to 
%                                    observed variables. NOTE: Delays are
%                                    reported as (real-valued) number of
%                                    time-steps.
%                    xDim         -- int; number of across-group latent 
%                                    variables
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
%     T       -- int; number of timesteps
%
% Outputs:
%
%     K_big   -- (xDim * numGroups * T) x (xDim * numGroups * T) array;
%                GP covariance matrix            
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.

xDim         = params.xDim;
yDims        = params.yDims;
numGroups    = length(yDims);
K_big        = zeros(xDim*numGroups*T);
K_big_inv    = zeros(xDim*numGroups*T);
mT           = numGroups*T;
logdet_K_big = 0;

Tdif = repmat(1:T,numGroups,1); % (m x T)
Tdif = repmat(Tdif(:)',mT,1) - repmat(Tdif(:),1,mT); % (mT x mT)
for i = 1:xDim
    Delayall = params.DelayMatrix(:,i); % (m x 1)
    Delaydif = repmat(Delayall,T,1);    % (mT x 1)
    Delaydif = repmat(Delaydif',mT,1) - repmat(Delaydif,1,mT); % (mT x mT)
    deltaT = Tdif - Delaydif; 
    deltaTsq = deltaT.^2;
    switch(params.covType)
        case 'rbf'
            temp = exp(-0.5*params.gamma(i)*deltaTsq);          
    end
    K_i = (1-params.eps(i))*temp + params.eps(i)*eye(mT); % (mT x mT)
    
    idx = i:xDim:xDim*numGroups*T;
    K_big(idx,idx) = K_i;
%     [K_big_inv(idx, idx), logdet_K] = invToeplitz(K_i);
%     
%     logdet_K_big = logdet_K_big + logdet_K;
end
