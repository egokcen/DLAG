function [K_big] = make_K_big_dlag(params,T)
%
% [K_big] = make_K_big_dlag(params, T)
%
% Description: Constructs full GP covariance matrix across all state 
%              dimensions (both within- and across-group) and timesteps.
%
% Arguments:
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
%                    if covType == 'sg'
%                        nu_across -- (1 x xDim_across) array; center
%                                     frequencies for spectral Gaussians;
%                                     convert to 1/time via 
%                                     nu_across./binWidth 
%                        nu_within -- (1 x numGroups) cell array; 
%                                     center frequencies for each group
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
%     T       -- int; number of timesteps
%
% Outputs:
%
%     K_big   -- (xDim * numGroups * T) x (xDim * numGroups * T) array;
%                GP covariance matrix            
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.
%     17 Apr 2020 -- Added 0-within-group dimension functionality
%     07 Jan 2021 -- Fixed error with indexing into K_big that only arose
%                    for numGroups > 2.
%     18 Feb 2023 -- Added spectral Gaussian compatibility.

% Initialize relevant variables
xDim_across  = params.xDim_across;
xDim_within  = params.xDim_within;
yDims        = params.yDims;
numGroups    = length(yDims);
K_big        = zeros((xDim_across*numGroups + sum(xDim_within))*T);

% Construct K for across-group latents
params_across.xDim = xDim_across;
params_across.yDims = yDims;
params_across.DelayMatrix = params.DelayMatrix;
params_across.gamma = params.gamma_across;
params_across.eps = params.eps_across;
params_across.covType = params.covType;
switch params.covType
    case 'sg'
        params_across.nu = params.nu_across;     
end
K_big_across = make_K_big_plusDelays(params_across,T);

% Construct K for within-group latents
K_big_within = cell(1,numGroups);
params_within.covType = params.covType;
for groupIdx = 1:numGroups
    % Leave K_big_within empty wherever xDim_within is 0
    if xDim_within(groupIdx) > 0
        params_within.xDim = xDim_within(groupIdx);
        params_within.gamma = params.gamma_within{groupIdx};
        params_within.eps = params.eps_within{groupIdx};
        params_within.C = params.C;
        switch params.covType
            case 'sg'
                params_within.nu = params.nu_within{groupIdx};     
        end
        % make_K_big is the same function as used for GPFA
        [K_big_within{groupIdx},~,~] = make_K_big(params_within,T,'xDim',xDim_within(groupIdx));
    end
end

% Fill in K_big with across- and within-group components
for t1 = 1:T
    % Time blocks
    bigBaseIdx1 = (xDim_across*numGroups + sum(xDim_within))*(t1-1) + 1;
    acrossBaseIdx1 = (xDim_across*numGroups)*(t1-1) + 1;
    for t2 = 1:T
        % Time blocks
        bigBaseIdx2 = (xDim_across*numGroups + sum(xDim_within))*(t2-1) + 1;
        acrossBaseIdx2 = (xDim_across*numGroups)*(t2-1) + 1;
        for i1 = 1:numGroups
            % Group blocks
            acrossIdx1 = acrossBaseIdx1 + (i1-1) * xDim_across;
            bigIdx1 = bigBaseIdx1 + (i1-1) * xDim_across + sum(xDim_within(1:i1-1));
            for i2 = 1:numGroups
                % Group blocks
                acrossIdx2 = acrossBaseIdx2 + (i2-1) * xDim_across;
                bigIdx2 = bigBaseIdx2 + (i2-1) * xDim_across + sum(xDim_within(1:i2-1));
                % Fill in across-group entries
                K_big(bigIdx1:bigIdx1+xDim_across-1, bigIdx2:bigIdx2+xDim_across-1) ...
                    = K_big_across(acrossIdx1:acrossIdx1+xDim_across-1,acrossIdx2:acrossIdx2+xDim_across-1);
            
                % Fill in within-group entries
                if i1 == i2
                    % Skip filling in K_big whenever xDim_within is 0
                    if xDim_within(i1) > 0
                        bigWithinIdx1 = bigIdx1+xDim_across:bigIdx1+xDim_across+xDim_within(i1)-1;
                        bigWithinIdx2 = bigIdx2+xDim_across:bigIdx2+xDim_across+xDim_within(i2)-1;
                        K_big(bigWithinIdx1,bigWithinIdx2) ...
                            = K_big_within{i1}((t1-1)*xDim_within(i1)+1:(t1-1)*xDim_within(i1)+xDim_within(i1), ...
                                               (t2-1)*xDim_within(i1)+1:(t2-1)*xDim_within(i1)+xDim_within(i1));
                    end
                end
            end
        end
    end
end


