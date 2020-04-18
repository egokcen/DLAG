function res = learnObsParams_dlag(seq, params, varargin)
%
% learnObsParams_dlag
%
% Description: Update parameters of the DLAG observation model given 
%              inferred latent states.
%
% Arguments:
%
%     Required:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
% 
%     params -- Structure containing DLAG model parameters. 
%               Contains the fields
% 
%               covType -- string; type of GP covariance (e.g., 'rbf')
%               gamma_across -- (1 x xDim_across) array; GP timescales
%                               in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
%               eps_across   -- (1 x xDim_across) GP noise variances
%               gamma_within -- (1 x numGroups) cell array; 
%                               GP timescales for each group
%               eps_within   -- (1 x numGroups) cell array;
%                               GP noise variances for each group
%               d            -- (yDim x 1) array; observation mean
%               C            -- (yDim x (numGroups*xDim)) array;
%                               mapping between low- and high-d spaces
%               R            -- (yDim x yDim) array; observation noise
%                               covariance 
%               DelayMatrix  -- (numGroups x xDim_across) array;
%                               delays from across-group latents to 
%                               observed variables. NOTE: Delays are
%                               reported as (real-valued) number of
%                               time-steps.
%               xDim_across  -- int; number of across-group latent 
%                               variables
%               xDim_within  -- (1 x numGroups) array; number 
%                               within-group latents in each group
%               yDims        -- (1 x numGroups) array; 
%                               dimensionalities of each observed group
%     
%     Optional:
%    
%     minVarFrac -- float; Set private variance floor, for entries of
%                   observation noise covariance matrix. (default: 0.01)
%
% Outputs:
%
%     res -- Structure containing updated observation model parameters:
%            C, d, and R (see 'params' above)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.   

    % Optional arguments
    minVarFrac = 0.01;
    assignopts(who, varargin);

    % Initialize other relevant variables
    xDim_across  = params.xDim_across;
    xDim_within  = params.xDim_within;
    yDims        = params.yDims;   % Number of features in each group
    yDim         = sum(yDims);     % Total number of features, all groups
    numGroups    = length(yDims);  % Number of groups
    N            = length(seq(:)); % Number of trials
    T            = [seq.T];        % List of trial lengths
    varFloor     = minVarFrac * diag(cov([seq.y]')); % Convert variance fraction to actual variance
    
    % Useful for extracting correct-sized blocks from matrices later
    block_idxs = get_block_idxs(yDims);

    % Solve for C and d together
    C = zeros(yDim,xDim_across*numGroups + sum(xDim_within));
    d = zeros(yDim,1);
    % This section is a bit tricky, since we have to deal with a variable
    % number of within-group dimensions for each group.
    for j = 1:numGroups
        currGroup = block_idxs{j};
        obs_idxs = currGroup(1):currGroup(2); % Current observation group
        % This line selects the correct within- and across-group latent
        % indexes for the current group, so we can access them in inferred
        % states and their posterior covariances.
        lat_idxs = (1+(j-1)*xDim_across+sum(xDim_within(1:j-1))):(j*xDim_across + sum(xDim_within(1:j))); % Current latent group
        sum_Pauto = zeros(xDim_across + xDim_within(j));
        for n = 1:N
            sum_Pauto = sum_Pauto + ...
                sum(seq(n).Vsm(lat_idxs,lat_idxs,:), 3) ...
                    + seq(n).xsm(lat_idxs,:) * seq(n).xsm(lat_idxs,:)';
        end
        Y = [seq.y]; Y = Y(obs_idxs,:); % Only take the current group
        Xsm = [seq.xsm]; Xsm = Xsm(lat_idxs,:);
        sum_yxtrans = Y * Xsm';
        sum_xall = sum(Xsm, 2);
        sum_yall = sum(Y, 2);
        
        term = [sum_Pauto sum_xall; sum_xall' sum(T)]; % (xDim+1) x (xDim+1)
        Cd   = ([sum_yxtrans sum_yall]) / term;        % yDim x (xDim+1)

        C(obs_idxs,lat_idxs) = Cd(:,1:length(lat_idxs));
        d(obs_idxs,1) = Cd(:,end);
    end
    % Collect outputs
    res.C = C;
    res.d = d;
    
    % Solve for R
    Y           = [seq.y];
    Xsm         = [seq.xsm];
    sum_yxtrans = Y * Xsm';
    sum_xall    = sum(Xsm, 2);
    sum_yall    = sum(Y, 2);
    
    if params.notes.RforceDiagonal
        % Constrain R to be diagonal
        sum_yytrans = sum(Y .* Y, 2);
        yd          = sum_yall .* res.d;
        term        = sum((sum_yxtrans - res.d * sum_xall') .* res.C, 2);
        r           = res.d.^2 + (sum_yytrans - 2*yd - term) / sum(T);
        
        % Set minimum private variance
        r           = max(varFloor, r);
        res.R = diag(r);    
    else
        % R is allowed to be a full covariance matrix
        sum_yytrans = Y * Y';
        yd          = sum_yall * res.d';
        term        = (sum_yxtrans - res.d * sum_xall') * res.C';
        R           = res.d * res.d' +...
                      (sum_yytrans - yd - yd' - term) / sum(T);
        
        res.R = (R+R')/2; % Ensure symmetry
    end

    if any(diag(res.R) == varFloor)
        fprintf('Warning: Private variance floor used for one or more observed dimensions in DLAG.\n');
    end
    
end 

