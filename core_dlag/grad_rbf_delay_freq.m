function [f,df] = grad_rbf_delay_freq(p,precomp,const)
%
% [f, df] = grad_rbf_delay_freq(p, precomp, const)  
%
% Description: Gradient computation for squared exponential GP 
%              hyperparameter (plus time delay) optimization, using
%              a frequency domain approximation. This function is called 
%              by minimize.m.
%
% Arguments:
%
%     p          -- variable with respect to which optimization is performed,
%                   where p = [ log(1 / timescale ^2), D(i,2),....D(i,M)]'
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
%     24 Jul 2023 -- Initial full revision.   
 

params = precomp.params;
numGroups = length(params.yDims);
df = zeros(size(p));
f = 0;
gamma = exp(p(1));
% Betaall = [0; p(2:end)];
% Delayall = params.maxDelay.*tanh(Betaall./2);
Delayall = [0; p(2:end)];
for j = 1:length(precomp.Tu)
    T = precomp.Tu(j).T;

    % Handle even and odd sequence lengths
    freqs = ((-floor(T/2):floor((T-1)/2))./T).';
      
    % Construct prior spectrum and inverse
    sqexp = (1 - const.eps).*exp(-0.5*((2*pi.*freqs).^2)./gamma);
    S = sqrt(2*pi/gamma) .* sqexp + const.eps;     % (T x 1) vector of diagonal elements
    S_inv = 1./S;
    logdet_S = sum(log(S));
            
    dS_dgamma = sqrt(pi/2) ...
        .* ((2*pi.*freqs).^2 .* gamma^(-5/2) - gamma^(-3/2)) ...
        .* sqexp;   % (xDim x 1) vector of diagonal elements
    Sinv_dSdgamma = S_inv .* dS_dgamma;
    XX_Sinv = precomp.Tu(j).XX .* S_inv.';
    df(1) = df(1) - 0.5.*trace((precomp.Tu(j).numTrials.*eye(T) - XX_Sinv).*Sinv_dSdgamma.');
    
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
