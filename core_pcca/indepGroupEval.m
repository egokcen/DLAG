function [LL, R2, MSE] = indepGroupEval(Ys, params, rGroups)
%
% [LL, R2, MSE] = indepGroupEval(Ys, params, rGroups)
%
% Description: Evaluate log likelihood of a pCCA model assuming all groups
%              are independent.
%
% Arguments:
%     Ys -- (1 x numGroups) cell array; list of data matrices 
%           {(y1Dim x N), (y2Dim x N), ...}
%     params.Rs -- (1 x numGroups) cell array; Blocks of uniqueness 
%                  matrix {(y1Dim x y1Dim), (y2Dim x y2Dim), ...}
%     params.ds -- (1 x numGroups) cell array; List of data means 
%                  {(y1Dim x 1), (y2Dim x 1), ...}
%     rGroups   -- (1 x 2) array; Each element specifies a group to be 
%                  included in the regression.
%
% Outputs:
%     LL  -- float; log likelihood of the data
%     R2  -- (1 x 2) array; R^2 in each pairwise direction 
%     MSE -- (1 x 2) array; mean-squared error in each pairwise direction
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Mar 2019 -- Initial full revision.

    numGroups = length(Ys);
    N = size(Ys{1},2);
    
    Rs = params.Rs;
    ds = params.ds;
    d = cat(1, params.ds{:});
    
    % Construct joint data matrix
    Y = cat(1,Ys{:});
    [yDim, N] = size(Y);
    
    Y0 = Y - repmat(d,1,N);
    cY = (1.0 / N) .* (Y0 * Y0');

    const = (-yDim / 2.0) * log(2.0 * pi);
    
    iRs = cell(1,numGroups);
    for groupIdx = 1:numGroups
        iRs{groupIdx} = inv(Rs{groupIdx});
    end
    iR = blkdiag(iRs{:}); % Inverse of R, (yDim x yDim) array

    % Compute the log likelihood
    ldM = sum(log(diag(chol(iR)))); % (1/2) log-determinant of iR
    LL = N*const + N*ldM - 0.5*N*sum(sum(iR .* cY));
    
    % Compute R^2 and MSE for the groups in rGroups
    numRGroups = length(rGroups);
    R2 = nan(1,numRGroups);
    MSE = nan(1,numRGroups);
    for rIdx = 1:numRGroups
        currSourceGroup = rGroups(rIdx);
        currTargetGroup = rGroups(setdiff(1:numRGroups,rIdx));
        Ytarget = Ys{currTargetGroup};
        % For pairwise regression, independent groups is
        % the same as predicting the mean every time
        d = ds{currTargetGroup};
        Ypred = repmat( d, [1 size(Ytarget,2)] );
        
        RSS = sum( sum( ( Ytarget - Ypred ).^2, 1 ) );
        TSS = sum( sum( ( Ytarget - repmat( mean(Ytarget,2), [1 size(Ytarget,2)] ) ).^2, 1 ) );
        R2(rIdx) = 1 - RSS / TSS;
        MSE(rIdx) = immse(Ypred, Ytarget);
    end
    
end