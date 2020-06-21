function [seq, LL, MSE, MSEorth, R2, R2orth] = denoise_dlag(seq, params)
%
% [seq, MSE, R2] = denoise_dlag(seq, params)
%
% Description: Denoise the data by inferring latent trajectories, and 
%              then projecting these trajectories back into observation 
%              space.
%
% Arguments:
%
%     seq     -- data structure, whose nth entry (corresponding to
%                the nth trial) has fields
%                    trialId      -- unique trial identifier
%                    T (1 x 1)    -- number of timesteps
%                    y (yDim x T) -- neural data
%
%     params  -- Structure containing DLAG model parameters at which EM
%                algorithm is initialized. Contains the fields
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
%                    xDim_within  -- (1 x numGroups) array; number of
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group%
% OUTPUTS:
%
%     seq     -- data structure with new fields yDenoised.
%                seq(n).yDenoised has the same dimensions as seq(n).y.
%     LL      -- float; data log-likelihood using the full DLAG model
%     MSE     -- float; mean-squared reconstruction error
%     MSEorth -- (xDim x 1) array; mean-squared error for the 
%                orthonormalized DLAG model.
%     R2      -- float; R^2 reconstruction error
%     R2orth  -- (xDim x 1) array; R^2 reconstruction error for the
%                orthonormalized DLAG model.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     11 Apr 2020 -- Initial full revision. 
%     08 Jun 2020 -- Updated to include denoising via the orthonormalized
%                    model.

    yDims = params.yDims;
    yDim = sum(yDims);
    numGroups = length(yDims);
    xDim = sum([numGroups*params.xDim_across, params.xDim_within]); % ALL latent dimensions

    % Initialize output structures
    for n = 1:length(seq)
        for m = 1:xDim 
            fn = sprintf('yDenoisedOrth%02d', m);
            seq(n).(fn) = nan(yDim, seq(n).T);
        end
    end
    
    MSEorth = nan(xDim,1);
    R2orth = nan(xDim,1);

    % Infer latent trajectories
    [seq, LL] = exactInferenceWithLL_dlag(seq, params, 'getLL', true);
    
    % Orthonormalize latents
    Xtemp = [seq.xsm];
    [Xorth, Corth] = orthogonalize(Xtemp, params.C);
    seq = segmentByTrial(seq, Xorth, 'xorth');
    
    % Project latents back into observation space
    for n = 1:length(seq)
        for m = 1:xDim
            fn = sprintf('yDenoisedOrth%02d', m);
            seq(n).(fn) = Corth(:,1:m) * seq(n).xorth(1:m,:) + params.d;
        end
    end
    
    % Now compute MSE and R^2
    Ytrue = [seq.y];
    % Orthonormalized DLAG performance
    for m = 1:xDim
        fn = sprintf('yDenoisedOrth%02d', m);
        Ypred = [seq.(fn)];
        % MSE
        MSEorth(m) = immse(Ypred, Ytrue);
        % R2
        RSS = sum( sum( ( Ytrue - Ypred ).^2, 1 ) );
        TSS = sum( sum( ( Ytrue - repmat( mean(Ytrue,2), [1 size(Ytrue,2)] ) ).^2, 1 ) );
        R2orth(m) = 1 - RSS / TSS;
    end
    % Full model performance
    MSE = MSEorth(end);
    R2 = R2orth(end);

end