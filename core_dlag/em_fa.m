function [estParams, LL] = em_fa(Y, xDim, varargin)
%
% [estParams, LL] = em_fa(Y, xDim, ...) 
%
% Description: Factor analysis and probabilistic PCA.
%
% Arguments: 
%
%     Required:
%
%     Y    -- (yDim x N) array; Observed data
%     xDim -- int; number of factors
%
%     Optional:
%
%     typ         -- string; 'fa' or 'ppca' (default: 'fa')
%     tolLL       -- float; stopping criterion for EM (default: 1e-8)
%     maxIters    -- int; maximum number of EM iterations (default: 1e8)
%     minVarFrac  -- float; fraction of overall data variance for each 
%                    observed dimension to set as the private variance 
%                    floor. (default: 0.01)
%     verbose     -- logical; that specifies whether to display status 
%                    messages (default: false)
%     parallelize -- logical; Here, this setting just determines which 
%                    status messages will be printed. (default: false)
%
% Outputs:
%
%     estParams.C  -- (yDim x xDim) array; factor loadings
%     estParams.R  -- (yDim x 1) array; diagonal of uniqueness matrix
%     estParams.d  -- (yDim x 1) array; data mean
%     LL           -- (1 x numIters) array; log likelihood at each EM 
%                     iteration
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%     Code adapted from fastfa.m by Byron Yu.
%
% Revision history:
%     27 Jun 2020 -- Initial full revision.

    typ         = 'fa';
    tolLL       = 1e-8; 
    maxIters    = 1e8;
    minVarFrac  = 0.01;
    verbose     = false;
    parallelize = false;
    assignopts(who, varargin);

    [yDim, N] = size(Y);

    % Initialization of parameters
    cY    = cov(Y', 1);
    if rank(cY) == yDim
        scale = exp(2*sum(log(diag(chol(cY))))/yDim);
    else
        % cX may not be full rank because N < xDim
        fprintf('WARNING in fastfa.m: Data matrix is not full rank.\n');
        r     = rank(cY);
        e     = sort(eig(cY), 'descend');
        scale = geomean(e(1:r));
    end
    C     = randn(yDim,xDim)*sqrt(scale/xDim);
    R    = diag(cY);
    d     = mean(Y, 2);

    varFloor = minVarFrac * diag(cY);  

    I     = eye(xDim);
    const = -xDim/2*log(2*pi);
    LLi   = 0; 
    LL    = [];
  
    for i = 1:maxIters
        % =======
        % E-step
        % =======
        iR  = diag(1./R);
        iRC = iR * C;  
        MM   = iR - iRC / (I + C' * iRC) * iRC';
        beta = C' * MM; % (xDim x yDim) array

        cY_beta = cY * beta'; % (yDim x xDim) array
        Ezz     = I - beta * C + beta * cY_beta;

        % Compute log likelihood
        LLold = LLi;    
        ldM   = sum(log(diag(chol(MM))));
        LLi   = N*const + N*ldM - 0.5*N*sum(sum(MM .* cY)); 
        if verbose && ~parallelize
          fprintf('EM iteration %3d of %d        lik %f\r', i, maxIters, LLi);
        end
        LL = [LL LLi];    

        % =======
        % M-step
        % =======
        C  = cY_beta / Ezz;
        R = diag(cY) - sum(cY_beta .* C, 2);

        if isequal(typ, 'ppca')
          R = mean(R) * ones(yDim, 1);
        end
        if isequal(typ, 'fa')
          % Set minimum private variance
          R = max(varFloor, R);
        end

        if i<=2
          LLbase = LLi;
        elseif (LLi < LLold)
          disp('VIOLATION');
        elseif ((LLi-LLbase) < (1+tolLL)*(LLold-LLbase))
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

    if any(R == varFloor)
        fprintf('Warning: Private variance floor used for one or more observed dimensions in FA.\n');
    end

    estParams.C = C;
    estParams.R = R;
    estParams.d = d;
