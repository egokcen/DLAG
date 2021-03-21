function snr = computeSNR_dlag(params)
%
% snr = computeSNR_dlag(params)
%
% Description: Compute the signal-to-noise ratio (SNR) for each group.
%              SNR is defined as the ratio of shared variance to private
%              noise variance (tr(CC')/tr(R)).
%
% Arguments:
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
%                    xDim_within  -- (1 x numGroups) array; number of 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
% Outputs:
%
%     snr -- (1 x numGroups) array; Signal-to-noise ratios for each group.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Mar 2021 -- Initial full revision.

% Initialize output structure
numGroups = length(params.yDims);
snr = nan(1,numGroups);

% Get observation parameters for each group
groupParams = partitionParams_dlag(params);

% Compute signal-to-noise ratios
for groupIdx = 1:numGroups
    R = groupParams{groupIdx}.R;
    C = groupParams{groupIdx}.C;
    snr(groupIdx) = trace(C*C') / trace(R);
end