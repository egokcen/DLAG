function [R2, MSE] = pred_dlag(seq, params)
%
% [R2, MSE] = pred_dlag(seq, params)
%
% Description: Performs leave-group-out prediction using an existing DLAG
%              model.
%
% Arguments:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
%
%     params  -- Structure containing DLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
%                    eps_across   -- (1 x xDim_across) GP noise variances
%                    gamma_within -- (1 x numGroups) cell array; 
%                                    GP timescales for each group
%                    eps_within   -- (1 x numGroups) cell array;
%                                    GP noise variances for each group
%                    d            -- (yDim x 1) array; observation mean
%                    C            -- (yDim x (numGroups*xDim)) array;
%                                    mapping between low- and high-d spaces
%                    R            -- (yDim x yDim) array; observation noise
%                                    covariance 
%                    DelayMatrix  -- (numGroups x xDim_across) array;
%                                    delays from across-group latents to 
%                                    observed variables. NOTE: Delays are
%                                    reported as (real-valued) number of
%                                    time-steps.
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
% Outputs:
%
%     MSE     -- float; leave-group-out mean-squared error
%     R2      -- float; leave-group-out R^2
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Jan 2023 -- Initial full revision.


yDims = params.yDims;
yDim = sum(yDims);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
numGroups = length(yDims);

obs_block_idxs = get_block_idxs(yDims); % Index rows of C
lat_block_idxs = get_block_idxs(xDim_across + xDim_within); % Index cols of C
xDim_all = numGroups*xDim_across + sum(xDim_within);

Ytrue = [seq.y]; % Stack all sequences

% Perform leave-group-out prediction
for groupIdx = 1:numGroups  
    
    % Indices of the current target group (responses)
    targetGroup = groupIdx;
    currTargetBlock = obs_block_idxs{targetGroup};
    targetIdxs = currTargetBlock(1):currTargetBlock(2);
    
    % Indices for delayed latents corresponding to target group
    currTargetLatBlock = lat_block_idxs{targetGroup};
    targetLatIdxs_all = currTargetLatBlock(1):currTargetLatBlock(2);
    targetLatIdxs = targetLatIdxs_all(1:xDim_across);
    
    % Indices of the current source group (predictors)
    sourceIdxs = setdiff(1:yDim,targetIdxs);
    
    % Keep only across-group latents
    
    for n = 1:length(seq)
        seqReg(n).trialId = seq(n).trialId;
        seqReg(n).T = seq(n).T;
        seqReg(n).y = seq(n).y(sourceIdxs,:);
    end
    
    if xDim_across > 0
        % Infer latents normally.
        paramsReg   = params;
        paramsReg.C = params.C(sourceIdxs,:);
        paramsReg.d = params.d(sourceIdxs);
        paramsReg.R = params.R(sourceIdxs,sourceIdxs);

        % Infer latents given only the source group
        seqReg = exactInferenceWithLL_dlag(seqReg, paramsReg, 'getLL', false);

        % Compute predicted activity of the target group
        for n = 1:length(seq)
            seq(n).ypred(targetIdxs,:) = params.C(targetIdxs,targetLatIdxs) * seqReg(n).xsm(targetLatIdxs,:) + params.d(targetIdxs);
        end
    else
        % Don't infer latents if xDim_across is 0.
        for n = 1:length(seq)
            seq(n).ypred(targetIdxs,:) = repmat(params.d(targetIdxs), 1, seq(n).T);
        end
    end
  
end

% Compute performance metrics
Ypred = [seq.ypred];
% MSE
MSE = immse(Ypred, Ytrue);
% R2
RSS = sum( sum( ( Ytrue - Ypred ).^2, 1 ) );
TSS = sum( sum( ( Ytrue - repmat( mean(Ytrue,2), [1 size(Ytrue,2)] ) ).^2, 1 ) );
R2 = 1 - RSS / TSS;
