function [estParams, LL] = em_pcca(Ys, xDim, varargin)
%
% [estParams, LL] = em_pcca(Ys, xDim, ...)
%
% Description: Fit a probabilistic canonical correlation analysis (pCCA)
%              model using the EM algorithm.
%
% Arguments:
%
%     Required:
%
%     Ys       -- (1 x numGroups) cell array; list of data matrices 
%                 {(y1Dim x N), (y2Dim x N), ...}
%     xDim     -- int; number of factors
%
%     Optional:
%
%     tolLL    -- float; stopping criterion for EM (default: 1e-8)
%     maxIters -- int; maximum number of EM iterations (default: 1e8)
%     verbose  -- boolean; specifies whether to display status messages
%                 (default: false)
%     parallelize -- logical; Here, this setting just determines which 
%                    status messages will be printed. (default: false)
%
% Outputs:
%     estParams.Cs -- (1 x numGroups) cell array; List of factor loadings 
%                     {(y1Dim x xDim), (y2Dim x xDim), ...}
%     estParams.Rs -- (1 x numGroups) cell array; Blocks of uniqueness 
%                     matrix {(y1Dim x y1Dim), (y2Dim x y2Dim), ...}
%     estParams.ds -- (1 x numGroups) cell array; List of data means 
%                     {(y1Dim x 1), (y2Dim x 1), ...}
%     LL           -- (1 x numIters) array; log likelihood at each EM 
%                     iteration
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Mar 2019 -- Initial full revision.
%     24 May 2020 -- Updated verbose output messages.

    tolLL       = 1e-8; 
    maxIters    = 1e8;
    verbose     = false;
    parallelize = false;
    assignopts(who, varargin);
    
    numGroups = length(Ys);
    % Get individual group sizes
    yDims = zeros(1,numGroups);
    for groupIdx = 1:numGroups
        yDims(groupIdx) = size(Ys{groupIdx},1); 
    end
    % Useful for extracting correct-sized blocks from matrices later
    block_idxs = get_block_idxs(yDims);
    block_mask = create_block_mask(yDims); % (yDim x yDim) array
    
    % Construct joint data matrix
    Y = cat(1,Ys{:});
    [yDim, N] = size(Y);
    
    % Initialize parameters
    cY = cov(Y', 1);
    if rank(cY) == yDim
        scale = exp(2*sum(log(diag(chol(cY))))/yDim); % Scale by determinant of cY
    else
        % cX may not be full rank because N < xDim
        fprintf('WARNING in pcca.m: Joint data matrix is not full rank.\n');
        r     = rank(cY);
        e     = sort(eig(cY), 'descend');
        scale = geomean(e(1:r));
    end
    C = randn(yDim,xDim)*sqrt(scale/xDim);
    Rs = cell(1,numGroups);
    for groupIdx = 1:numGroups
        Rs{groupIdx} = cov(Ys{groupIdx}',1);
    end
    R = blkdiag(Rs{:});
    d = mean(Y,2);
    
    % Constants
    I = eye(xDim);
    const = (-yDim / 2.0) * log(2.0 * pi);
    LLi = 0.0;
    LL = [];
    
    for i = 1:maxIters
        % =======
        % E-step
        % =======
        iRs = cell(1,numGroups);
        for groupIdx = 1:numGroups
            iRs{groupIdx} = inv(Rs{groupIdx});
        end
        iR = blkdiag(iRs{:}); % Inverse of R, (yDim x yDim) array
        iR = 0.5 * (iR + iR'); % Ensure symmetry
        iRC = iR * C;
        MM = iR - iRC / (I + C' * iRC) * iRC'; % Inverse of (CC' + R), (yDim x yDim array)
        beta = C' * MM; % (xDim x yDim) array
        
        cY_beta = cY * beta'; % (yDim x xDim) array
        Exx = I - beta * C + beta * cY_beta;
        
        % Compute the log likelihood
        LLold = LLi;
        ldM = sum(log(diag(chol(MM)))); % (1/2) log-determinant of MM
        LLi = N*const + N*ldM - 0.5*N*sum(sum(MM .* cY));
        if verbose && ~parallelize
            fprintf('EM iteration %3d of %d        lik %f\r', i, maxIters, LLi);
        end
        LL = [LL LLi]; 
        
        % =======
        % M-step
        % =======
        % Intermediate computations
        C = cY_beta / Exx;
        R = cY - cY_beta * C'; % Full update
        % Enforce R to be symmetric (helps with numerical stability too)
        R = 0.5 * (R + R');
        R = R .* block_mask;
        Rs = mat2blocks(R, block_idxs);
        
        % Check stopping conditions or errors
        if i <= 2
            LLbase = LLi;
        elseif (LLi < LLold)
            disp('VIOLATION');
        elseif ((LLi - LLbase) < (1+tolLL)*(LLold - LLbase))
            break;
        end
         
    end
  
    if verbose && ~parallelize
        if length(LL) < maxIters
            fprintf('LL converged after %d EM iterations.\n', length(LL));
        else
            fprintf('Fitting stopped after maxIters (%d) was reached.\n', maxIters);
        end 
    end
    
    % Convert joint params into lists of group params
    estParams.Cs = cell(1,numGroups);
    estParams.ds = cell(1,numGroups);
    for groupIdx = 1:numGroups
        currGroup = block_idxs{groupIdx};
        estParams.Cs{groupIdx} = C(currGroup(1):currGroup(2),:);
        estParams.ds{groupIdx} = d(currGroup(1):currGroup(2));
    end
    estParams.Rs = Rs;
    
end