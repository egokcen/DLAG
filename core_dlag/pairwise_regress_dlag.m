function [seq, MSE, MSEorth, R2, R2orth] = pairwise_regress_dlag(seq, params, rGroups)
%
% [seq, MSE, MSEorth, R2, R2orth] = pairwise_regress_dlag(seq, params, rGroups)
%
% Description: Performs reduced-rank regression between two groups using
%              the DLAG model.
%              NOTE: Within-group (target) latent activity is excluded.
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
%     rGroups -- (1 x 2) array; Each element specifies a group to be 
%                included in the regression.
%
% Outputs:
%
%     seq     -- test data structure with new fields yregOrthXX, where XX are
%                the number of predictive dimensions used.  
%                seq(n).yregOrthXX has the same dimensions as seq(n).y.
%     MSE     -- structure with the following fields:
%                indiv -- (1 x 2) array; mean-squared error in each 
%                         pairwise direction
%                joint -- float; mean-squared error across both groups
%     MSEorth -- structure with the following fields:
%                indiv -- (xDim_across+1 x 2) array; mean-squared error in 
%                         each direction for orthonoralized DLAG predictions
%                joint -- (xDim_across+1 x 1) array; mean-squared error 
%                         across both groups for orthonormalized DLAG 
%                         predictions
%     R2      -- structure with the following fields:
%                indiv -- (1 x 2) array; R^2 in each pairwise direction
%                joint -- float; R^2 computed across both groups jointly
%     R2orth  -- structure with the following fields:
%                indiv -- (xDim_across+1 x 2) array; R^2 error in each 
%                         direction for orthonormalized DLAG predictions
%                joint -- (xDim_across+1 x 1) array; R^2 computed across 
%                         both groups jointly for orthonormalized DLAG 
%                         predictions.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.
%     08 Jun 2020 -- Updated to include joint performance metrics.
%     27 Jun 2020 -- Added 0-across-group dimension functionality.
%     20 Oct 2021 -- Overhauled to use predictiveProjection_dlag.m


yDims = params.yDims;
yDim = sum(yDims);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;

obs_block_idxs = get_block_idxs(yDims); % Index rows of C
lat_block_idxs = get_block_idxs(xDim_across + xDim_within); % Index cols of C

numRGroups = length(rGroups); % Really just 2, for pairwise regression

% Initialize output structures
% Predicted sequences
for n = 1:length(seq)
    for m = 0:xDim_across
        fn = sprintf('yregOrth%02d', m);
        seq(n).(fn) = nan(yDim, seq(n).T);
    end
end

Ytrue_joint = [seq.y]; % Stack all sequences

% Performance metrics
MSE.indiv = nan(1,numRGroups);
MSEorth.indiv = nan(xDim_across+1,numRGroups);
MSEorth.joint = nan(xDim_across+1,1);
R2.indiv = nan(1,numRGroups);
R2orth.indiv = nan(xDim_across+1,numRGroups);
R2orth.joint = nan(xDim_across+1,1);
  
% Compute predictions 
for i = 1:numRGroups  
    
    % Indices of the current source group (predictors)
    currSourceGroup = rGroups(i);
    currSourceBlock = obs_block_idxs{currSourceGroup};
    sourceIdxs = currSourceBlock(1):currSourceBlock(2);
    
    % Indices for delayed latents corresponding to source group
    currSourceLatBlock = lat_block_idxs{currSourceGroup};
    sourceLatIdxs = currSourceLatBlock(1):currSourceLatBlock(2);
    
    % Indices of the current target group (responses)
    currTargetGroup = rGroups(setdiff(1:numRGroups,i));
    currTargetBlock = obs_block_idxs{currTargetGroup};
    targetIdxs = currTargetBlock(1):currTargetBlock(2);
    
    % Indices for delayed latents corresponding to target group
    currTargetLatBlock = lat_block_idxs{currTargetGroup};
    targetLatIdxs = currTargetLatBlock(1):currTargetLatBlock(2);
    targetLatIdxs = targetLatIdxs(1:xDim_across); % Keep only across-group latents
    
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
        
        % Project latents onto orthonormal basis for predictive subspace
        % (in target group)
        seqReg = predictiveProjection_dlag(seqReg, params, ...
            'groupIdxs', [currSourceGroup currTargetGroup], ...
            'orth', true);
        
        % Compute predictive modes
        [~, ~, ~, V, ~] = predictiveModes_dlag(params, ...
            'groupIdxs',[currSourceGroup currTargetGroup]);


        % Compute predicted activity of the target group
        for n = 1:length(seq)
            for m = 0:xDim_across
                fn = sprintf('yregOrth%02d', m);
                if m > 0
                    seq(n).(fn)(targetIdxs,:) = V{2}(:,1:m) * seqReg(n).xpred(xDim_across+(1:m),:) + params.d(targetIdxs);
                else
                    seq(n).(fn)(targetIdxs,:) = repmat(params.d(targetIdxs), 1, seq(n).T); 
                end
            end
        end
    else
        % Don't infer latents if xDim_across is 0.
        for n = 1:length(seq)
            fn = sprintf('yregOrth%02d', xDim_across);
            seq(n).(fn)(targetIdxs,:) = repmat(params.d(targetIdxs), 1, seq(n).T);
        end
    end
        
    % Compute individual performance metrics
    Ytrue_indiv = Ytrue_joint(targetIdxs,:);
    % Orthonormalized DLAG performance
    for m = 0:xDim_across
        fn = sprintf('yregOrth%02d', m);
        Ypred = [seq.(fn)];
        Ypred = Ypred(targetIdxs,:);
        % MSE
        MSEorth.indiv(m+1,i) = immse(Ypred, Ytrue_indiv);
        % R2
        RSS = sum( sum( ( Ytrue_indiv - Ypred ).^2, 1 ) );
        TSS = sum( sum( ( Ytrue_indiv - repmat( mean(Ytrue_indiv,2), [1 size(Ytrue_indiv,2)] ) ).^2, 1 ) );
        R2orth.indiv(m+1,i) = 1 - RSS / TSS;
    end
    % Full model performance
    MSE.indiv(i) = MSEorth.indiv(end,i);
    R2.indiv(i) = R2orth.indiv(end,i);
  
end

% Compute joint performance metrics
% Orthonormalized model performance
for m = 0:xDim_across
    fn = sprintf('yregOrth%02d', m);
    Ypred = [seq.(fn)];
    % MSE
    MSEorth.joint(m+1,1) = immse(Ypred, Ytrue_joint);
    % R2
    RSS = sum( sum( ( Ytrue_joint - Ypred ).^2, 1 ) );
    TSS = sum( sum( ( Ytrue_joint - repmat( mean(Ytrue_joint,2), [1 size(Ytrue_joint,2)] ) ).^2, 1 ) );
    R2orth.joint(m+1,1) = 1 - RSS / TSS;
end
% Full model performance
MSE.joint = MSEorth.joint(end);
R2.joint = R2orth.joint(end);
