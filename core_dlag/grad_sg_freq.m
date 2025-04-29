function [f, df] = grad_sg_freq(p, precomp, const)
%
% [f, df] = grad_sg_freq(p, precomp, const)  
%
% Description: Gradient computation for spectral Gaussian GP hyperparameter
%              optimization, using a frequency domain approximation. This 
%              function is called by minimize.m.
%
% Arguments:
%
%     p       -- variable with respect to which optimization is performed,
%                where p = [ log(1 / timescale ^2), center freq.]'
%     precomp -- structure containing precomputations
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
nu = p(2);
dE_dgamma = 0;
dE_dnu = 0;
f = 0;
df = zeros(size(p));
for j = 1:length(precomp.Tu)
    T = precomp.Tu(j).T;

    % Handle even and odd sequence lengths
    freqs = ((-floor(T/2):floor((T-1)/2))./T).';
    freqsmin = freqs-nu;
    freqspls = freqs+nu;
    sqfreqsmin = (2*pi.*(freqsmin)).^2;
    sqfreqspls = (2*pi.*(freqspls)).^2;

    % Construct prior spectrum and inverse
    sqexpmin = exp(-0.5*sqfreqsmin./gamma);
    sqexppls = exp(-0.5*sqfreqspls./gamma);
    S = (1-const.eps).*sqrt(pi/(2*gamma)) .* (sqexpmin + sqexppls) + const.eps; % (xDim x 1) vector of diagonal elements
    S_inv = 1./S;
    logdet_S = sum(log(S));
            
    % gamma
    dS_dgamma = (1-const.eps) * sqrt(pi/8).*( ...
        -gamma^(-3/2) .* (sqexpmin + sqexppls) ...
        + gamma^(-5/2) .* (sqfreqsmin.*sqexpmin + sqfreqspls.*sqexppls));   % (xDim x 1) vector of diagonal elements

    Sinv_dSdgamma = S_inv .* dS_dgamma;
    XX_Sinv = precomp.Tu(j).XX .* S_inv;
    dE_dgamma = dE_dgamma - 0.5 * sum((precomp.Tu(j).numTrials - XX_Sinv) .* Sinv_dSdgamma);
        
    % nu
    dS_dnu = (1-const.eps) * sqrt(8*pi^5/gamma^3) ...
        .* (freqsmin.*sqexpmin - freqspls.*sqexppls);   % (xDim x 1) vector of diagonal elements
    Sinv_dSdnu = S_inv .* dS_dnu;
    dE_dnu = dE_dnu - 0.5 * sum((precomp.Tu(j).numTrials - XX_Sinv) .* Sinv_dSdnu);

    f = f - 0.5 * precomp.Tu(j).numTrials * logdet_S - 0.5 * sum(XX_Sinv);
end

f  = -f;
df(1) = -dE_dgamma * gamma; % df/d(log(gamma))
df(2) = -dE_dnu;