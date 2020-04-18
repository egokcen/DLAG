function [seq, MSE, R2] = denoise_dlag(seq, params)
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
%     MSE     -- float; mean-squared reconstruction error
%     R2      -- float; R^2 reconstruction error
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     11 Apr 2020 -- Initial full revision. 

    yDims = params.yDims;
    yDim = sum(yDims);

    fn = sprintf('yDenoised');
    for n = 1:length(seq)
      seq(n).(fn) = nan(yDim, seq(n).T);
    end

    % Infer latent trajectories
    seq = exactInferenceWithLL_dlag(seq, params, 'getLL', false);
    
    % Project latents back into observation space
    for n = 1:length(seq)
        seq(n).(fn) = params.C * seq(n).xsm + params.d;
    end
    
    % Now compute MSE and R^2
    Ytrue = [seq.y];
    Ypred = [seq.(fn)];
    % MSE
    MSE = immse(Ypred, Ytrue);
    % R2
    RSS = sum( sum( ( Ytrue - Ypred ).^2, 1 ) );
    TSS = sum( sum( ( Ytrue - repmat( mean(Ytrue,2), [1 size(Ytrue,2)] ) ).^2, 1 ) );
    R2 = 1 - RSS / TSS;

end