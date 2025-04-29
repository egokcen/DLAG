function k = transformed_gpcov(V, params, binWidth, varargin)
%
% k = transformed_gpcov(V, params, binWidth, ...)
%
% Description: Given a coupled set of basis vectors, V, compute GP
%              covariance functions for projections that lie along those
%              directions.
%
% Arguments:
%
%     Required:
%
%     V       -- (1 x numGroups) cell array; V{i} is a (yDims(i) x r)
%                array, where r is the number of basis vectors.
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
%     normalize -- logical; true to compute correlation functions (i.e.,
%                    normalized); false to compute covariance functions
%                    (default: true)
%
% Outputs:
%
%     k -- (1 x r) cell array; k{r} is a (numGroups x numGroups) cell
%          array, and k{r}{i,j} is the cross-covariance function between
%          group i and group j. k{r}{i,j} is an anonymous function.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Feb 2022 -- Initial full revision.
%     21 Feb 2023 -- Added spectral Gaussian compatibility. Transitioned to
%                    anonymous functions. Ignoring private variance now.

showPlot = true;
stepres = 0.1;
maxlag = 4*binWidth;
normalize = true;
assignopts(who,varargin);

numGroups = length(V);
r = size(V{1},2);
xDim_within = params.xDim_within;
xDim_across = params.xDim_across;
lagSteps = -maxlag:stepres:maxlag;

% Initialize output structure
k = cell(1,r);
for j = 1:r
    k{j} = cell(numGroups); 
end

% Get observation parameters for each group
groupParams = partitionParams_dlag(params);
Ca = cell(1,numGroups);
Cw = cell(1,numGroups);
R = cell(1,numGroups);
for groupIdx = 1:numGroups
    R{groupIdx} = groupParams{groupIdx}.R;
    acrossParams = getSubsetParams_dlag(groupParams{groupIdx}, 1:xDim_across, cell(1,numGroups));
    Ca{groupIdx} = acrossParams.C; 
    withinParams = getSubsetParams_dlag(groupParams{groupIdx},[], {1:xDim_within(groupIdx)});
    Cw{groupIdx} = withinParams.C;
end

% Compute prior covariance functions
[ka_prior, kw_prior] = gpcov(params, binWidth, ...
                             'normalize', true, ...
                             'showPlot', false);

% Compute transformed covariance functions
for j = 1:r
    for groupIdx1 = 1:numGroups
        for groupIdx2 = 1:numGroups
            k{j}{groupIdx1,groupIdx2} = @(t) 0;
            for xIdx = 1:xDim_across
                k{j}{groupIdx1,groupIdx2} = @(t) k{j}{groupIdx1,groupIdx2}(t) ...
                    + ka_prior{groupIdx1,groupIdx2}{xIdx}(t) ...
                    .*(V{groupIdx1}(:,j)'*Ca{groupIdx1}(:,xIdx)*Ca{groupIdx2}(:,xIdx)'*V{groupIdx2}(:,j));
            end
            if groupIdx1 == groupIdx2
                for xIdx = 1:xDim_within(groupIdx1)
                    k{j}{groupIdx1,groupIdx2} = @(t) k{j}{groupIdx1,groupIdx2}(t) ...
                    + kw_prior{groupIdx1}{xIdx}(t) ...
                    .*(V{groupIdx1}(:,j)'*Cw{groupIdx1}(:,xIdx)*Cw{groupIdx1}(:,xIdx)'*V{groupIdx1}(:,j)); 
                end
            end
        end
    end
end

if normalize
    % Normalize covariance functions to get correlation functions
    k_norm = k;
    for j = 1:r
        for groupIdx1 = 1:numGroups
            for groupIdx2 = 1:numGroups
                k_norm{j}{groupIdx1,groupIdx2} ...
                    = @(t) k{j}{groupIdx1,groupIdx1}(0).^(-0.5) ...
                      .* k{j}{groupIdx1,groupIdx2}(t) ...
                      .* k{j}{groupIdx2,groupIdx2}(0).^(-0.5);
            end
        end
    end
    k = k_norm;
end

if showPlot
    if normalize
        ylbl = 'Correlation';
    else
        ylbl = 'Covariance';
    end
    
    for j = 1:r
        figure;
        for groupIdx1 = 1:numGroups
            for groupIdx2 = 1:numGroups
                k_curr = k{j}{groupIdx1,groupIdx2}(lagSteps);
                subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                hold on;
                [maxval, maxdelay] = max(abs(k_curr));
                if k_curr(maxdelay) < 0
                    k_curr = -1.*k_curr; 
                end
                minval = min(k_curr);
                minval = 0;
                line([lagSteps(maxdelay) lagSteps(maxdelay)], [minval 1.05*max([1 maxval])], 'color', 'r', 'linestyle', '--');
                line([0 0], [minval 1.05*max([1 maxval])], 'color', 'k', 'linestyle', '--');
                plot(lagSteps, k_curr, 'k-', 'linewidth', 1.5);
                axis square;
                xlabel('Time lag (ms)');
                ylabel(ylbl);
                axis([min(lagSteps) max(lagSteps) minval 1.05*max([1 maxval])]);
            end
        end
    end
end