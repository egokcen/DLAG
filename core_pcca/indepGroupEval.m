function [LL, R2, MSE] = indepGroupEval(Ys, params, rGroups)
%
% [LL, R2, MSE] = indepGroupEval(Ys, params, rGroups)
%
% Description: Evaluate log likelihood of a pCCA model assuming all groups
%              are independent.
%
% Arguments:
%
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
%
%     LL  -- float; log likelihood of the data
%     R2  -- structure with the following fields:
%            indiv -- (1 x 2) array; R^2 in each pairwise direction
%            joint -- float; R^2 computed across both groups jointly
%                        (leave-group-out prediction)
%     MSE -- structure with the following fields:
%            indiv -- (1 x 2) array; mean-squared error in each 
%                      pairwise direction
%            joint -- float; mean-squared error across both groups
%                     (leave-group-out prediction)
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Mar 2019 -- Initial full revision.
%     25 Feb 2022 -- Added leave-group-out (joint) prediction metrics.

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
    Ypred = cell(1,numRGroups);
    Ytrue = cell(1,numRGroups);
    R2.indiv = nan(1,numRGroups);
    MSE.indiv = nan(1,numRGroups);
    for rIdx = 1:numRGroups
        currSourceGroup = rGroups(rIdx);
        currTargetGroup = rGroups(setdiff(1:numRGroups,rIdx));
        Ytarget = Ys{currTargetGroup};
        Ytrue{rIdx} = Ytarget;
        % For pairwise regression, independent groups is
        % the same as predicting the mean every time
        d = ds{currTargetGroup};
        Ypred{rIdx} = repmat( d, [1 size(Ytarget,2)] );
        
        RSS = sum( sum( ( Ytarget - Ypred{rIdx} ).^2, 1 ) );
        TSS = sum( sum( ( Ytarget - repmat( mean(Ytarget,2), [1 size(Ytarget,2)] ) ).^2, 1 ) );
        R2.indiv(rIdx) = 1 - RSS / TSS;
        MSE.indiv(rIdx) = immse(Ypred{rIdx}, Ytarget);
    end
    
    % Compute joint performance metrics
    Ytrue_joint = cat(1,Ytrue{:});
    Ypred_joint = cat(1,Ypred{:});

    MSE.joint = immse(Ypred_joint, Ytrue_joint);
    RSS = sum( sum( ( Ytrue_joint - Ypred_joint ).^2, 1 ) );
    TSS = sum( sum( ( Ytrue_joint - repmat( mean(Ytrue_joint,2), [1 size(Ytrue_joint,2)] ) ).^2, 1 ) );
    R2.joint = 1 - RSS / TSS;
    
end