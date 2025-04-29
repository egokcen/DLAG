function [eap_across, eap_within] = gppowers(params, binWidth)
%
% [eap_across, eap_within] = gppowers(params, binWidth)
%
% Description: Compute the theoretical (cross) covariances (i.e., expected
%              average power) for each latent via their (cross) power
%              spectral densities.
%
% Arguments:
%
%     Required:
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
%     binWidth   -- float; resolution (sample period or bin width), in
%                   units of time, at which params were estimated.
%
% Outputs:
%
%     eap_across -- (numGroups x numGroups x xDim) array; 
%                   eap_across(j,i,k) is the expected average power for 
%                   across-area latent j, between group i and group k
%     eap_within -- (1 x numGroups) cell array; eap_within{i} is a 
%                   (1 x xDim_within(i)) array with the expected average
%                   power of each within-area latent.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     13 Feb 2023 -- Initial full revision.

numGroups = length(params.yDims);
xDim_within = params.xDim_within;
xDim_across = params.xDim_across;

% Convert GP params to same units as binWidth
gp_params = getGPparams_dlag(params,binWidth);

% Compute expected average powers (covariances)
eap_across = nan(numGroups,numGroups,xDim_across); % Across-area amplitude functions
eap_within = cell(1,numGroups);           % Within-area amplitude functions
for groupIdx1 = 1:numGroups
    for groupIdx2 = 1:numGroups
        for xIdx = 1:xDim_across
            D = (gp_params.DelayMatrix(groupIdx2,xIdx) - gp_params.DelayMatrix(groupIdx1,xIdx));
            ka = @(f) sqrt(2*pi*(gp_params.tau_across(xIdx))^2) ...
                .*exp(-0.5*(gp_params.tau_across(xIdx)^2).*(2*pi.*f).^2) ...
                .* exp(1i.*D .* (2*pi.*f));
            eap_across(groupIdx1,groupIdx2,xIdx) = integral(ka,-Inf,Inf);
        end
        if groupIdx1 == groupIdx2
            eap_within{groupIdx1} = nan(1,xDim_within(groupIdx1));
            for xIdx = 1:xDim_within(groupIdx1)
                kw = @(f) sqrt(2*pi*gp_params.tau_within{groupIdx1}(xIdx)^2) ...
                    .*exp(-0.5*(gp_params.tau_within{groupIdx1}(xIdx)^2).*(2*pi.*f).^2);
                eap_within{groupIdx1}(xIdx) = integral(kw,-Inf,Inf);
            end
        end
    end
end