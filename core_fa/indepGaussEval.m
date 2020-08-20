function [LL, sse] = indepGaussEval(Y, params)
%
% [LL, sse] = indepGaussEval(Y, params)
%
% Description: Evaluate log likelihood and sum-of-squared prediction error 
%              given an independent Gaussian model
%
% Arguments: 
%
%     Y           -- (yDim x N) array; data matrix
%     estParams.d -- (yDim x 1) array; mean of each observed dimension
%     estParams.R -- (yDim x 1) array; variance of each observed dimension
%
% Outputs:
%
%     LL  -- float; log likehood of data
%     sse -- float; sum-of-squared prediction error
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

    [yDim, N] = size(Y);

    R   = params.R;
    d   = params.d;

    Yc2  = bsxfun(@minus, Y, d).^2;
    Ystd = bsxfun(@rdivide, Yc2, R);

    LL = -0.5 * (N*yDim*log(2*pi) + N*sum(log(R)) + sum(Ystd(:)));

    sse = sum(Yc2(:));