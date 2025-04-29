function [f,df] = grad_delay_freq(p,precomp,const)
%
% [f, df] = grad_delay_freq(p, precomp, const)  
%
% Description: Gradient computation for GP delay optimization, using a 
%              frequency domain approximation. This function is called by 
%              minimize.m.
%
% Arguments:
%
%     p          -- variable with respect to which optimization is 
%                   performed, where p = [D(m,1); ....; D(m,p)]
%     precomp    -- structure containing precomputations
%     const      -- structure containing parameters that stay constant
%                   during this optimization
%
% Outputs:
%
%     f          -- value of portion of lower bound that depends on p
%     df         -- gradient of f at p    
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     31 Aug 2023 -- Initial full revision.   
 
xDim = length(p);
df = zeros(size(p));
f = 0;
D = const.maxDelay.*tanh(p./2); % Compute unconstrained gradients
for j = 1:length(precomp.Tu)
    T = precomp.Tu(j).T;
    
    % Handle even and odd sequence lengths
    freqs = ((-floor(T/2):floor((T-1)/2))./T).';
    ifreqs = 1i*2*pi*freqs;

    % Construct time delay operators
    Q = exp(-D*ifreqs.'); % (xDim x T)
    QQXXCRinvC = repmat(permute(conj(Q),[1 3 2]),[1 xDim 1]) ...
              .* repmat(permute(Q,[3 1 2]),[xDim 1 1]) ...
              .* precomp.Tu(j).XXCRinvC;

    QYXRinvC = repmat(permute(Q,[3 1 2]),[size(precomp.Tu(j).YXRinvC,1) 1 1]) ...
            .* precomp.Tu(j).YXRinvC;

    % Gradient update
    df = df + permute(...
                  real(sum(...
                      repmat(permute(ifreqs,[2 3 1]),1,xDim,1).*(...
                          sum(QQXXCRinvC,1) - sum(QYXRinvC,1)...
                      ),...
                  3)),...
              [2 1 3]);
    % Objective update
    f = f - real(0.5.*sum(QQXXCRinvC,'all') - sum(QYXRinvC,'all'));

end
f = -f;
df = -df.*(const.maxDelay/2).*(sech(p./2)).^2; % Respect constraints
