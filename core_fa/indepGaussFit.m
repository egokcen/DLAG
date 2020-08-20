function estParams = indepGaussFit(Y)
%
% estParams = indepGaussFit(Y) 
%
% Description: Fit parameters of an independent Gaussian model.
%
% Arguments:
%
%     Y -- (yDim x N) array; data matrix
%
% Outputs:
%
%     estParams.d -- (yDim x 1) array; mean of each observed dimension
%     estParams.R -- (yDim x 1) array; variance of each observed dimension
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

    estParams.R = var(Y, 1, 2);
    estParams.d  = mean(Y, 2);
