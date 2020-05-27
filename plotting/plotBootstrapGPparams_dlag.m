function plotBootstrapGPparams_dlag(params,bootParams,binWidth,rGroups,varargin)
%
% plotBootstrapGPparams_dlag(params,bootParams,binWidth,rGroups,...)
%
% Description: Plot DLAG Delay Matrix and latent timescales, along with 
%              bootstrap confidence intervals (see bootstrapGPparams.m). 
%              Delays and across-group timescales are plotted as ordered 
%              pairs, only for the reference groups specified in rGroups.
%
% Arguments: 
%
%     Required:
%
%     params    -- Structure containing DLAG model parameters. 
%                  Contains the fields
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
%     bootParams  -- Structure containing bootstrapped DLAG GP parameters.
%                    Contains the fields
%
%                    tau_across.upper -- (1 x xDim_across) array; 
%                                         confidence interval upper bound
%                                         on across-group GP timescales  
%                    tau_across.lower -- (1 x xDim_across) array; 
%                                         confidence interval lower bound
%                                         on across-group GP timescales 
%                    tau_within.upper -- (1 x numGroups) cell array; 
%                                        confidence interval upper bound
%                                        on GP timescales for each group
%                    tau_within.lower -- (1 x numGroups) cell array; 
%                                        confidence interval lower bound
%                                        on GP timescales for each group
%                    DelayMatrix.upper -- (numGroups x xDim_across) array;
%                                        confidence interval upper bound on 
%                                        the delay matrix, converted to 
%                                        units of time.
%                    DelayMatrix.lower -- (numGroups x xDim_across) array;
%                                        confidence interval lower bound on 
%                                        the delay matrix, converted to 
%                                        units of time.
%
%     binWidth   -- float; bin width or sample period (in e.g., ms)
%     rGroups    -- (1 x 2) array; Indexes of groups to get relative
%                   delays. rGroups(1) is the reference group.
%
%     Optional:
%
%     plotAcross -- logical; set to true to plot across-group GP params
%                   (default: true)
%     plotWithin -- logical; set to true to plot within-group GP params
%                   (default: true)
%     units      -- string; units of time of binWidth (for labels)
%                   (default: '')
%
% Outputs:
%     
%     None. (But creates figures)
%              
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 May 2020 -- Initial full revision.

% Set optional arguments
plotAcross = true;
plotWithin = true;
units = '';
assignopts(who,varargin);

xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
yDims = params.yDims;
numGroups = length(yDims);
colors = generateColors(); % Generate custom plotting colors

% Convert GP params into units of time
gp_params = getGPparams_dlag(params, binWidth);

% Take only the delays associated with the groups in rGroups
delays = gp_params.DelayMatrix(rGroups(2),:) ...
       - gp_params.DelayMatrix(rGroups(1),:);
bootDelays.upper = bootParams.DelayMatrix.upper(rGroups(2),:) ...
                 - bootParams.DelayMatrix.upper(rGroups(1),:);
bootDelays.lower = bootParams.DelayMatrix.lower(rGroups(2),:) ...
                 - bootParams.DelayMatrix.lower(rGroups(1),:);

% Figure out plot limits
maxDelay = max([delays bootDelays.upper bootDelays.lower]);
minDelay = min([delays bootDelays.upper bootDelays.lower]);

all_tau = [gp_params.tau_across ...
           bootParams.tau_across.upper ...
           bootParams.tau_across.lower ...
           gp_params.tau_within{:} ...
           bootParams.tau_within.upper{:} ...
           bootParams.tau_within.lower{:}];

maxTau = max(all_tau);

% Determine number of subplots, based on optional arguments
numPlot = sum([plotAcross plotWithin.*(sum(xDim_within > 0))]);

if numPlot <= 0
    % Reach here if plotAcross is false and all xDim_within are 0
    fprintf('plotGPparams_dlag: No GP parameters to plot.\n')
else
    figure;
end
plotIdx = 0;

% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

% Across-group GP params
if plotAcross
    plotIdx = plotIdx + 1;
    subplot(1,numPlot,plotIdx);
    hold on;
    xlabel(sprintf('Delay from area %d to area %d%s',rGroups(1),rGroups(2),units));
    ylabel(sprintf('Across-group GP timescale%s', units));
    xlim([minDelay-10,maxDelay+10]);
    ylim([0,maxTau+10]);
    line([0 0], [0 maxTau+10], ...
         'Color', colors.grays{6}, ...
         'linestyle', '--', ...
         'linewidth', 1.5);
    errorbar(delays, gp_params.tau_across, ...
             gp_params.tau_across - bootParams.tau_across.lower, ...
             bootParams.tau_across.upper - gp_params.tau_across, ...
             delays - bootDelays.lower, ...
             bootDelays.upper - delays, ...
             'o', ...
             'color', colors.grays{6}, ...
             'linestyle', 'none', ...
             'linewidth', 1.5, ...
             'MarkerFaceColor', colors.reds{3});
    hold off;
end

% Within-group GP timescales
if plotWithin
    for groupIdx = 1:numGroups
        % Don't try to plot anything for groups with 0 within-group latents
        if xDim_within(groupIdx) > 0
            plotIdx = plotIdx + 1;
            subplot(1,numPlot,plotIdx);
            hold on;
            xlim([0,xDim_within(groupIdx)+1]); 
            ylim([0,maxTau+10]);
            errorbar([1:xDim_within(groupIdx)], gp_params.tau_within{groupIdx}, ...
                     gp_params.tau_within{groupIdx} - bootParams.tau_within.lower{groupIdx}, ...
                     bootParams.tau_within.upper{groupIdx} - gp_params.tau_within{groupIdx}, ...
                     'o', ...
                     'color', colors.grays{6}, ...
                     'linestyle', 'none', ...
                     'linewidth', 1.5, ...
                     'MarkerFaceColor', colors.reds{3});
            ylabel(sprintf('GP timescale%s', units));
            set(gca,'XTick',1:xDim_within(groupIdx));
            xlabel(sprintf('Within-group latents, group %d', groupIdx));
            hold off;
        end
    end
end