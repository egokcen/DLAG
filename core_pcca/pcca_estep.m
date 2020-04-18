function [X, LL] = pcca_estep(Ys, params)
%
% [X, LL] = pcca_estep(Ys, params)
%
% Description: Compute the low-dimensional points and data likelihoods
%              using a previously learned PCCA model.
%
% Arguments:
%     Ys -- (1 x numGroups) cell array; list of data matrices 
%           {(y1Dim x N), (y2Dim x N), ...}
%     params -- learned PCCA parameteres (structure with fields Cs, Rs, ds)
%   
% Outputs:
%     X.mean -- (xDim x N) array; posterior mean
%     X.cov  -- (xDim x xDim) array; posterior covariance
%     LL     -- log-likelihood of data
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Mar 2019 -- Initial full revision.

    numGroups = length(Ys);
    xDim = size(params.Cs{1}, 2);
    
    C = cat(1, params.Cs{:});
    Rs = params.Rs;
    R = blkdiag(Rs{:});
    d = cat(1, params.ds{:});
    
    % Construct joint data matrix
    Y = cat(1,Ys{:});
    [yDim, N] = size(Y);
    
    Y0 = Y - repmat(d,1,N);
    cY = (1.0 / N) .* (Y0 * Y0');
    
    % Constants
    I = eye(xDim);
    const = (-yDim / 2.0) * log(2.0 * pi);
    
    iRs = cell(1,numGroups);
    for groupIdx = 1:numGroups
        iRs{groupIdx} = inv(Rs{groupIdx});
    end
    iR = blkdiag(iRs{:}); % Inverse of R, (yDim x yDim) array
    iR = 0.5 * (iR + iR'); % Ensure symmetry
    iRC = iR * C;
    MM = iR - iRC / (I + C' * iRC) * iRC'; % Inverse of (CC' + R), (yDim x yDim array)
    beta = C' * MM; % (xDim x yDim) array

    X.mean = beta * Y0; % (xDim x N) array
    X.cov = I - beta * C; % (xDim x xDim) array; same for all observations

    % Compute the log likelihood
    ldM = sum(log(diag(chol(MM)))); % (1/2) log-determinant of MM
    LL = N*const + N*ldM - 0.5*N*sum(sum(MM .* cY));
    
end