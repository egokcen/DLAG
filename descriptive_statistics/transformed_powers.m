function eap = transformed_powers(V, params, binWidth, groupIdx)
%
% eap = transformed_powers(V, params, binWidth, groupIdx)
%
% Description: Given a coupled set of basis vectors, V, evaluate
%              theoretical (cross) covariances between modes from 
%              theoretical (cross) power spectral densities.
%
% Arguments:
%
%     Required:
%
%     V       -- (yDims(i) x r) array; r is the number of basis vectors.
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
%     groupIdx   -- int; evaluate a specific group given by groupIdx
%
% Outputs:
%
%     eap -- (r x r) array; eap(i,j) is the theoretical expected average
%            power shared between mode i and mode j.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     13 Feb 2023 -- Initial full revision.

r = size(V,2);
xDim_within = params.xDim_within(groupIdx);
xDim_across = params.xDim_across;
xDim = xDim_across + xDim_within;

% Initialize output structures
eap = nan(r);

% Convert GP params to same units as binWidth
gp_params = getGPparams_dlag(params,binWidth);

% Get observation parameters for the specified group
groupParams = partitionParams_dlag(params);
R = groupParams{groupIdx}.R;
acrossParams = getSubsetParams_dlag(groupParams{groupIdx}, 1:xDim_across, {[]});
Ca = acrossParams.C; 
withinParams = getSubsetParams_dlag(groupParams{groupIdx},[], {1:xDim_within});
Cw = withinParams.C;

% Compute prior spectral densities
spec_a_prior = cell(1,xDim_across);
for xIdx = 1:xDim_across
    D = (gp_params.DelayMatrix(groupIdx,xIdx) - gp_params.DelayMatrix(groupIdx,xIdx));
    spec_a_prior{xIdx} ...
        = @(f) sqrt(2*pi*gp_params.tau_across(xIdx)^2) ...
        .*exp(-0.5*(gp_params.tau_across(xIdx)^2).*(2*pi.*f).^2) ...
        .* exp(1i.*D .* (2*pi.*f));
end
spec_w_prior = cell(1,xDim_within);
for xIdx = 1:xDim_within
    spec_w_prior{xIdx} ...
        = @(f) sqrt(2*pi*gp_params.tau_within{groupIdx}(xIdx)^2) ...
        .*exp(-0.5*(gp_params.tau_within{groupIdx}(xIdx)^2).*(2*pi.*f).^2);
end

% Compute expected average powers
for r1 = 1:r
    for r2 = 1:r
        spec_mixed = @(f) 0;
        for xIdx = 1:xDim_across
            spec_mixed = @(f) spec_mixed(f) + (V(:,r1)'*Ca(:,xIdx))*(Ca(:,xIdx)'*V(:,r2)) .* spec_a_prior{xIdx}(f);
        end

        for xIdx = 1:xDim_within
            spec_mixed = @(f) spec_mixed(f) + (V(:,r1)'*Cw(:,xIdx))*(Cw(:,xIdx)'*V(:,r2)) .* spec_w_prior{xIdx}(f);
        end
        
        eap(r1,r2) = integral(spec_mixed,-Inf,Inf);
    end
end