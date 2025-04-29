function [f, df] = grad_sg(p, precomp, const)
%
% [f, df] = grad_sg(p, precomp, const)  
%
% Description: Gradient computation for parameters of the spectral Gaussian
%              GP kernel. This function is called by minimize.m.
%
% Arguments:
%
%     p          -- variable with respect to which optimization is performed,
%                   where p = [ log(1 / timescale ^2), center freq.]'
%     precomp    -- structure containing precomputations
% 
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
%     19 Feb 2023 -- Initial full revision.
  
  Tall = precomp.Tall;
  Tmax = max(Tall);
  
  temp1        = (1-const.eps) * exp(-exp(p(1)) / 2 * precomp.difSq); % Tmax x Tmax
  temp2        = temp1 .* cos(2*pi*p(2) * precomp.dif);                 
  Kmax         = temp2 + const.eps * eye(Tmax);
  dKdgamma_max = -0.5 * temp2 .* precomp.difSq;
  dKdnu_max    = -2*pi*precomp.dif .* temp1 .* sin(2*pi*p(2) * precomp.dif);
  
  dEdgamma = 0;
  dEdnu    = 0;
  f        = 0;
  df       = zeros(size(p));
  for j = 1:length(precomp.Tu)
    T     = precomp.Tu(j).T;
    Thalf = ceil(T/2);
    mkr = ceil(0.5 * T^2);
            
    [Kinv, logdet_K] = invToeplitz(Kmax(1:T, 1:T));
    
    % Gamma
    KinvM_g     = Kinv(1:Thalf,:) * dKdgamma_max(1:T,1:T);  % Thalf x T
    KinvMKinv_g = (KinvM_g * Kinv)';                         % Thalf x T
    
    dg_KinvM_g  = diag(KinvM_g);
    tr_KinvM_g  = 2 * sum(dg_KinvM_g) - rem(T, 2) * dg_KinvM_g(end);
    
    dEdgamma = dEdgamma - 0.5 * precomp.Tu(j).numTrials * tr_KinvM_g...
    + 0.5 * precomp.Tu(j).PautoSUM(1:mkr) * KinvMKinv_g(1:mkr)'... 
    + 0.5 * precomp.Tu(j).PautoSUM(end:-1:mkr+1) * KinvMKinv_g(1:(T^2-mkr))';

    % nu
    KinvM_n     = Kinv(1:Thalf,:) * dKdnu_max(1:T,1:T);  % Thalf x T
    KinvMKinv_n = (KinvM_n * Kinv)';                     % Thalf x T
    
    dn_KinvM_n  = diag(KinvM_n);
    tr_KinvM_n  = 2 * sum(dn_KinvM_n) - rem(T, 2) * dn_KinvM_n(end);
    
    dEdnu = dEdnu - 0.5 * precomp.Tu(j).numTrials * tr_KinvM_n...
    + 0.5 * precomp.Tu(j).PautoSUM(1:mkr) * KinvMKinv_n(1:mkr)'... 
    + 0.5 * precomp.Tu(j).PautoSUM(end:-1:mkr+1) * KinvMKinv_n(1:(T^2-mkr))';
        
    f = f - 0.5 * precomp.Tu(j).numTrials * logdet_K...
    - 0.5 * precomp.Tu(j).PautoSUM(:)' * Kinv(:);
  end
  
  f  = -f;
  % exp(p) is needed because we're computing gradients with
  % respect to log(gamma), rather than gamma
  df(1) = -dEdgamma * exp(p(1));
  df(2) = -dEdnu;
