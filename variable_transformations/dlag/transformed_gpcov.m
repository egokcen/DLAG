function k = transformed_gpcov(V, params, binWidth, maxlag, varargin)
%
% k = transformed_gpcov(V, params, binWidth, maxlag, ...)
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
%     maxlag     -- float; maximum time lag to consider when computing
%                   covariance function; same units as binWidth.
%
%     Optional:
%
%     showPlot -- logical; set true to plot all transformed covariance
%                 functions (default: true)
%     stepres  -- float; resolution of the computed covariance 
%                 functions, in the same units of time as binWidth
%                 (default: 0.1).
%     computeCorr -- logical; true to compute correlation functions (i.e.,
%                    normalized); false to compute covariance functions
%                    (default: true)
%
% Outputs:
%
%     k -- (1 x r) cell array; k{r} is a (numGroups x numGroups) cell
%          array, and k{r}{i,j} is the cross-covariance function between
%          group i and group j.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Feb 2022 -- Initial full revision.

showPlot = true;
stepres = 0.1;
computeCorr = true;
assignopts(who,varargin);

numGroups = length(V);
r = size(V{1},2);
xDim_within = params.xDim_within;
xDim_across = params.xDim_across;
lagSteps = -maxlag:stepres:maxlag;
T = length(lagSteps);

k_cov = cell(1,r);

% Convert GP params to same units as binWidth
gp_params = getGPparams_dlag(params,binWidth);

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
ka_prior = cell(numGroups,numGroups);
kw_prior = cell(1,numGroups);
for groupIdx1 = 1:numGroups
    for groupIdx2 = 1:numGroups
        ka_prior{groupIdx1,groupIdx2} = cell(1,xDim_across);
        for xIdx = 1:xDim_across
            D = (gp_params.DelayMatrix(groupIdx2,xIdx) - gp_params.DelayMatrix(groupIdx1,xIdx));
            ka_prior{groupIdx1,groupIdx2}{xIdx} ...
                = (1 - params.eps_across(xIdx)).*exp(-0.5*(gp_params.tau_across(xIdx)^(-2)).*(lagSteps-D).^2);
            ka_prior{groupIdx1,groupIdx2}{xIdx}((lagSteps - D) == 0) ...
                = ka_prior{groupIdx1,groupIdx2}{xIdx}((lagSteps - D) == 0) ...
                  + params.eps_across(xIdx);
        end
        if groupIdx1 == groupIdx2
            kw_prior{groupIdx1} = cell(1,xDim_within(groupIdx1));
            for xIdx = 1:xDim_within(groupIdx1)
                kw_prior{groupIdx1}{xIdx} ...
                    = (1 - params.eps_within{groupIdx1}(xIdx)).*exp(-0.5*(gp_params.tau_within{groupIdx1}(xIdx)^(-2)).*(lagSteps).^2); 
                kw_prior{groupIdx1}{xIdx}(lagSteps == 0) ...
                    = kw_prior{groupIdx1}{xIdx}(lagSteps == 0) ...
                      + params.eps_within{groupIdx1}(xIdx);
            end
        end
    end
end

% Compute transformed covariance functions
for j = 1:r
    for groupIdx1 = 1:numGroups
        for groupIdx2 = 1:numGroups
            k_cov{j}{groupIdx1,groupIdx2} = zeros(1,T);
            for xIdx = 1:xDim_across
                k_cov{j}{groupIdx1,groupIdx2} = k_cov{j}{groupIdx1,groupIdx2} ...
                    + ka_prior{groupIdx1,groupIdx2}{xIdx} ...
                    .*(V{groupIdx1}(:,j)'*Ca{groupIdx1}(:,xIdx)*Ca{groupIdx2}(:,xIdx)'*V{groupIdx2}(:,j));
            end
            if groupIdx1 == groupIdx2
                for xIdx = 1:xDim_within(groupIdx1)
                    k_cov{j}{groupIdx1,groupIdx2} = k_cov{j}{groupIdx1,groupIdx2} ...
                    + kw_prior{groupIdx1}{xIdx} ...
                    .*(V{groupIdx1}(:,j)'*Cw{groupIdx1}(:,xIdx)*Cw{groupIdx1}(:,xIdx)'*V{groupIdx1}(:,j)); 
                end
                k_cov{j}{groupIdx1,groupIdx2}(lagSteps == 0) = k_cov{j}{groupIdx1,groupIdx2}(lagSteps == 0) ...
                    + V{groupIdx1}(:,j)'*R{groupIdx1}*V{groupIdx1}(:,j);
            end
            if max(k_cov{j}{groupIdx1,groupIdx2}) < 0
                k_cov{j}{groupIdx1,groupIdx2} = -1.*k_cov{j}{groupIdx1,groupIdx2};
            end
        end
    end
end

if computeCorr
    % Normalize covariance functions to get correlation functions
    k = k_cov;
    for j = 1:r
        for groupIdx1 = 1:numGroups
            for groupIdx2 = 1:numGroups
                k{j}{groupIdx1,groupIdx2} ...
                    = (k_cov{j}{groupIdx1,groupIdx1}(lagSteps == 0).^(-0.5)) ...
                      .* k_cov{j}{groupIdx1,groupIdx2} ...
                      .* (k_cov{j}{groupIdx2,groupIdx2}(lagSteps == 0).^(-0.5));
            end
        end
    end
else
    k = k_cov; 
end

if showPlot
    for j = 1:r
        figure;
        for groupIdx1 = 1:numGroups
            for groupIdx2 = 1:numGroups
                k_curr = k{j}{groupIdx1,groupIdx2};
                subplot(numGroups,numGroups,(groupIdx1-1)*numGroups+groupIdx2);
                hold on;
                plot(lagSteps, k_curr, 'k-', 'linewidth', 1.5);
                [maxval, maxdelay] = max(k_curr);
                line([lagSteps(maxdelay) lagSteps(maxdelay)], [0 1.05*max([1 maxval])], 'color', 'r', 'linestyle', '--');
                line([0 0], [0 1.05*max([1 maxval])], 'color', 'k', 'linestyle', '--');
                axis square;
                xlabel('Time lag (ms)');
                ylabel('Correlation');
                axis([min(lagSteps) max(lagSteps) 0 1.05*max([1 maxval])]);
            end
        end
    end
end