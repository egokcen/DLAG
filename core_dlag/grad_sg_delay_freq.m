function [f,df] = grad_sg_delay_freq(p,precomp,const)
%
% [f, df] = grad_sg_delay_freq(p, precomp, const)  
%
% Description: Gradient computation for spectral Gaussian GP hyperparameter
%              optimization, using a frequency domain approximation. This 
%              function is called by minimize.m.
%
% Arguments:
%
%     p          -- variable with respect to which optimization is performed,
%                   where p = [ log(1 / timescale ^2), D(i,2),...,D(i,q),
%                   center freq]'
%     precomp    -- structure containing precomputations
%     const      -- structure containing parameters that stay constant
%                   during this optimization
%
% Outputs:
%
%     f          -- value of objective function E[log P({x},{y})] at p
%     df         -- gradient at p    
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     25 Jul 2023 -- Initial full revision.   
 

params = precomp.params;
numGroups = length(params.yDims);
df = zeros(size(p));
f = 0;
gamma = exp(p(1));
% Betaall = [0; p(2:numGroups)];
% Delayall = params.maxDelay.*tanh(Betaall./2);
Delayall = [0; p(2:numGroups)];
nu = p(end);
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
    XX_Sinv = precomp.Tu(j).XX .* S_inv.';
    df(1) = df(1) - 0.5.*trace((precomp.Tu(j).numTrials.*eye(T) - XX_Sinv).*Sinv_dSdgamma.');
        
    % nu
    dS_dnu = (1-const.eps) * sqrt(8*pi^5/gamma^3) ...
        .* (freqsmin.*sqexpmin - freqspls.*sqexppls);   % (xDim x 1) vector of diagonal elements
    Sinv_dSdnu = S_inv .* dS_dnu;
    df(end) = df(end) - 0.5.*trace((precomp.Tu(j).numTrials.*eye(T) - XX_Sinv).*Sinv_dSdnu.');
    
    % Delays
    % dDelayall_dBetaall = (params.maxDelay/2).*(sech(Betaall./2)).^2;
    ifreqs = 1i*2*pi*freqs;
    for k = 2:length(Delayall)%2:length(Betaall)
        QXYRC = exp(-ifreqs.*Delayall(k)).*precomp.Tu(j).XYRC(:,k); % (T x 1);
        dE_dDk = sum(real(-ifreqs.*QXYRC));
        df(k) = df(k) + dE_dDk;%*dDelayall_dBetaall(k);
        f = f + sum(real(QXYRC));
    end
      
    f = f - 0.5 * precomp.Tu(j).numTrials * logdet_S - 0.5 * trace(XX_Sinv);  
end
f = -f;
df(1) = df(1)*gamma; % df/d(log(gamma))
df = -df;
