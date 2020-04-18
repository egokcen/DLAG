function [f,df] = grad_betgam_delay(p,precomp,const)
%
% [f, df] = grad_betgam_delay(p, precomp, const)  
%
% Description: Gradient computation for GP timescale and delay 
%              optimization. This function is called by minimize.m.
%
% Arguments:
%
%     p          -- variable with respect to which optimization is performed,
%                   where p = [ log(1 / timescale ^2), Tau(i,2),....Tau(i,q)]'
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
%     18 Mar 2020 -- Initial full revision.   
 

  params = precomp.params;
  numGroups = length(params.yDims);
  df = zeros(size(p));
  f = 0;
  gamma = exp(p(1));
  Betaall = [0; p(2:end)];
  Delayall = 2*params.maxDelay./(1+exp(-Betaall)) ...
                - params.maxDelay;
  for j = 1:length(precomp.Tu)
      T = precomp.Tu(j).T;
      mT = numGroups*T;
      dK_dBetak = zeros(mT);
      
      Delaydif = repmat(Delayall,T,1); 
      Delaydif = repmat(Delaydif',mT,1) - repmat(Delaydif,1,mT);
      deltaT = (precomp.Tdif(1:mT,1:mT) - Delaydif); 
      deltaTsq = deltaT.^2;
      temp = (1-const.eps)*exp(-(gamma/2) * deltaTsq);
      dtemp = gamma * temp .* deltaT;
      
      K = temp + const.eps*eye(mT);
      KinvPautoSUM = K\precomp.Tu(j).PautoSUM;
      dE_dK = -0.5*(precomp.Tu(j).numTrials*eye(mT) - KinvPautoSUM)/K;
      dK_dgamma = -0.5*temp.*deltaTsq;
      dE_dgamma = dE_dK(:)' * dK_dgamma(:);
      df(1) = df(1) + dE_dgamma;
      
      exp_Beta = exp(-Betaall);
      dDelayall_dBetaall = -2*params.maxDelay*exp_Beta./((1+exp_Beta).^2);
      for k = 2:length(Betaall)
          idx = k:numGroups:numGroups*T;
          dK_dBetak(:,idx) = -dtemp(:,idx)*dDelayall_dBetaall(k);
          dK_dBetak(idx,:) = dtemp(idx,:)*dDelayall_dBetaall(k);
          dK_dBetak(idx,idx) = 0;
          dE_dBetak = dE_dK(:)' * dK_dBetak(:);
          df(k) = df(k) + dE_dBetak;
          dK_dBetak(:,idx) = 0;
          dK_dBetak(idx,:) = 0;
      end
      
      f = f - 0.5*precomp.Tu(j).numTrials*logdet(K) -0.5*trace(KinvPautoSUM); 
  end
f = -f;
df(1) = df(1)*gamma; % df/d(log(gamma))
df = -df;
