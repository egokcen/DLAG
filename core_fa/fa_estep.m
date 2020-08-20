function [X, LL] = fastfa_estep(Y, params)
%
% [X, LL] = fa_estep(Y, params)
%
% Description: Compute the low-dimensional points and data likelihoods 
%              using a previously learned FA or PPCA model.
%
% Arguments:
%
%     Y         -- (yDim x N) array; data matrix
%     params.C  -- (yDim x xDim) array; factor loadings
%     params.R  -- (yDim x 1) array; diagonal of uniqueness matrix
%     params.d  -- (yDim x 1) array; data mean
%
% Outputs:
%
%     X.mean -- (xDim x N) array; posterior mean
%     X.cov  -- (xDim x xDim) array; posterior covariance, which is the 
%               same for all data
%     LL     -- float; log-likelihood of data
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

  [yDim, N] = size(Y);
  xDim      = size(params.C, 2);
    
  C     = params.C;
  R     = params.R;
  d     = params.d;
  
  Yc   = bsxfun(@minus, Y, d);
  YcYc = Yc * Yc';

  I = eye(xDim);
    
  const = -yDim/2*log(2*pi);
    
  iR  = diag(1./R);
  iRC = iR * C;    
  MM   = iR - iRC / (I + C' * iRC) * iRC';
  beta = C' * MM; % xDim x yDim
    
  X.mean = beta * Yc; % xDim x N
  X.cov  = I - beta * C; % xDim x xDim; same for all observations

  LL = N*const + 0.5*N*logdet(MM) - 0.5 * sum(sum(MM .* YcYc));
