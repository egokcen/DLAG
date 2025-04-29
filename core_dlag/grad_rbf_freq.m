function [f, df] = grad_rbf_freq(p, precomp, const)
%
% [f, df] = grad_rbf_freq(p, precomp, const)  
%
% Description: Gradient computation for squared exponential GP 
%              hyperparameter optimization, using a frequency domain 
%              approximation. This function is called by minimize.m.
%
% Arguments:
%
%     p       -- variable with respect to which optimization is performed,
%                where p = log(1 / timescale ^2)
%     precomp -- structure containing precomputations
%     const   -- structure containing parameters that stay constant
%                during this optimization
%
% Outputs:
%
%     f       -- value of objective function E[log P({x},{y})] at p
%     df      -- gradient at p    
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     21 Jul 2023 -- Initial full revision. 
%     13 Aug 2024 -- Updated gradient and lower bound with much more
%                    efficient computation of XX.
    
gamma = exp(p(1));
dE_dgamma = 0;
f = 0;
for j = 1:length(precomp.Tu)
    T = precomp.Tu(j).T;

    % Handle even and odd sequence lengths
    freqs = ((-floor(T/2):floor((T-1)/2))./T).';

    % Construct prior spectrum and inverse
    sqexp = (1 - const.eps).*exp(-0.5*((2*pi.*freqs).^2)./gamma);
    S = sqrt(2*pi/gamma) .* sqexp + const.eps;     % (xDim x 1) vector of diagonal elements
    S_inv = 1./S;
    logdet_S = sum(log(S));
            
    dS_dgamma = sqrt(pi/2) ...
        .* ((2*pi.*freqs).^2 .* gamma^(-5/2) - gamma^(-3/2)) ...
        .* sqexp;   % (xDim x 1) vector of diagonal elements
    Sinv_dSdgamma = S_inv .* dS_dgamma;
    XX_Sinv = precomp.Tu(j).XX .* S_inv;
    dE_dgamma = dE_dgamma - 0.5 * sum((precomp.Tu(j).numTrials - XX_Sinv) .* Sinv_dSdgamma);
        
    f = f - 0.5 * precomp.Tu(j).numTrials * logdet_S - 0.5 * sum(XX_Sinv);
end

f  = -f;
df = -dE_dgamma * gamma; % df/d(log(gamma))
