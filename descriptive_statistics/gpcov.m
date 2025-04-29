function [ka, kw] = gpcov(params, binWidth, varargin)
%
% [ka, kw] = gpcov(params, binWidth, ...)
%
% Description: Compute and (optionally) plot GP covariance functions over
%              a range of time lags specified by maxlag.
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
%                    xDim_within  -- (1 x numGroups) array; number of 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
%     binWidth   -- float; resolution (sample period or bin width), in
%                   units of time, at which params were estimated.
%
%     Optional:
%
%     showPlot -- logical; set true to plot all transformed covariance
%                 functions (default: true)
%     maxlag   -- float; maximum time lag to consider when computing
%                 covariance function; same units as binWidth.
%                 (default: 4*binWidth)
%     stepres  -- float; resolution of the computed covariance 
%                 functions, in the same units of time as binWidth
%                 (default: 0.1).
%     normalize -- logical; set true to compute cross-correlation functions,
%                  else cross-covariance functions. (default: true)
%
% Outputs:
%
%     ka -- (numGroups x numGroups) cell array; ka{i,j} is a 
%           (1 x xDim_across) cell array, and ka{i,j}{r} is the cross-
%           covariance function between group i and group j, given by
%           latent r. ka{i,j}{r} in an anonymous function.
%     kw -- (1 x numGroups) cell array; kw{i} is a (1 x xDim_within{i}) 
%           cell array, and kw{i}{r} is the within-group covariance 
%           function for group i, given by latent r. kw{i}{r} in an 
%           anonymous function.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     30 Dec 2022 -- Initial full revision.
%     18 Feb 2023 -- Added spectral Gaussian compatibility.
%     20 Feb 2023 -- Transitioned to anonymous functions, and ignoring
%                    GP noise variance now.

showPlot = true;
stepres = 0.1;
maxlag = 4*binWidth;
normalize = true;
assignopts(who,varargin);

numGroups = length(params.yDims);
xDim_within = params.xDim_within;
xDim_across = params.xDim_across;
lagSteps = -maxlag:stepres:maxlag;

% Convert GP params to same units as binWidth
gp_params = getGPparams_dlag(params,binWidth);
groupParams = partitionParams_dlag(params);

% Compute prior covariance functions
ka = cell(numGroups,numGroups);
kw = cell(1,numGroups);
for groupIdx1 = 1:numGroups
    for groupIdx2 = 1:numGroups
        ka{groupIdx1,groupIdx2} = cell(1,xDim_across);
        for xIdx = 1:xDim_across
            D = (gp_params.DelayMatrix(groupIdx2,xIdx) - gp_params.DelayMatrix(groupIdx1,xIdx));
            switch params.covType
                case 'rbf'
                    ka{groupIdx1,groupIdx2}{xIdx} ...
                        = @(t) exp(-0.5*(gp_params.tau_across(xIdx)^(-2)).*(t-D).^2);
                case 'sg'
                    ka{groupIdx1,groupIdx2}{xIdx} ...
                        = @(t) exp(-0.5*(gp_params.tau_across(xIdx)^(-2)).*(t-D).^2) ...
                          .*cos(2*pi*gp_params.nu_across(xIdx).*(t-D));
            end
            if ~normalize
                ka{groupIdx1,groupIdx2}{xIdx} = @(t) norm(groupParams{groupIdx1}.C(:,xIdx)) ...
                    .* ka{groupIdx1,groupIdx2}{xIdx}(t) .* norm(groupParams{groupIdx2}.C(:,xIdx));
            end
        end
        if groupIdx1 == groupIdx2
            kw{groupIdx1} = cell(1,xDim_within(groupIdx1));
            for xIdx = 1:xDim_within(groupIdx1)
                switch params.covType
                    case 'rbf'
                        kw{groupIdx1}{xIdx} ...
                            = @(t) exp(-0.5*(gp_params.tau_within{groupIdx1}(xIdx)^(-2)).*(t).^2);
                    case 'sg'
                        kw{groupIdx1}{xIdx} ...
                            = @(t) exp(-0.5*(gp_params.tau_within{groupIdx1}(xIdx)^(-2)).*(t).^2) ...
                              .*cos(2*pi*gp_params.nu_within{groupIdx1}(xIdx).*(t));
                end
                if ~normalize
                    kw{groupIdx1}{xIdx} = @(t) norm(groupParams{groupIdx1}.C(:,xDim_across+xIdx)) ...
                        .*kw{groupIdx1}{xIdx}(t) .* norm(groupParams{groupIdx2}.C(:,xDim_across+xIdx));
                end
            end
        end
    end
end

if showPlot
    % Across-group covariance functions
    if xDim_across > 0
        for xIdx = 1:xDim_across
            figure;
            for groupIdx1 = 1:numGroups
                for groupIdx2 = 1:numGroups
                    k_curr = ka{groupIdx1,groupIdx2}{xIdx}(lagSteps);
                    subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                    hold on;
                    [maxval, maxdelay] = max(k_curr);
                    minval = min(k_curr);
                    minval = 0;
                    line([lagSteps(maxdelay) lagSteps(maxdelay)], [minval 1.05*max([1 maxval])], 'color', 'r', 'linestyle', '--');
                    line([0 0], [minval 1.05*max([1 maxval])], 'color', 'k', 'linestyle', '--');
                    plot(lagSteps, k_curr, 'k-', 'linewidth', 1.5);
                    axis square;
                    xlabel('Time lag (ms)');
                    ylabel('Correlation');
                    axis([min(lagSteps) max(lagSteps) minval 1.05*max([1 maxval])]);
                end
            end
        end
    end
    
    % Within-group covariance functions
    for groupIdx = 1:numGroups
        if xDim_within(groupIdx) > 0
            figure;
            for xIdx = 1:xDim_within(groupIdx)
                k_curr = kw{groupIdx}{xIdx}(lagSteps);
                subplot(1,xDim_within(groupIdx),xIdx);
                hold on;
                [maxval, maxdelay] = max(k_curr);
                minval = min(k_curr);
                minval = 0;
                line([lagSteps(maxdelay) lagSteps(maxdelay)], [0 1.05*max([1 maxval])], 'color', 'r', 'linestyle', '--');
                line([0 0], [minval 1.05*max([1 maxval])], 'color', 'k', 'linestyle', '--');
                plot(lagSteps, k_curr, 'k-', 'linewidth', 1.5);
                axis square;
                xlabel('Time lag (ms)');
                ylabel('Correlation');
                axis([min(lagSteps) max(lagSteps) minval 1.05*max([1 maxval])]);
            end
        end
    end
end