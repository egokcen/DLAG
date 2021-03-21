function gp_params = plotGPparams_dlag(params,binWidth,rGroups,varargin)
%
% gp_params = plotGPparams_dlag(params,binWidth,rGroups,...)
%
% Description: Plot DLAG Delay Matrix and latent timescales. Delays and 
%              across-group timescales are plotted as ordered pairs, only
%              for the reference groups specified in rGroups.
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
%     prom -- Shade plotted latents according to prominence. prom 
%             contains the following (relevant) fields:
%             across.LL.raw -- (1 x numDimGroups) array; prominence of
%                              dimGroups_across(i) evaluated on raw data, 
%                              measured by decrease in log-likelihood 
%                              relative to the full model.
%             across.VE.raw -- (1 x numDimGroups) array;prominence of
%                              dimGroups_across(i) evaluated on raw data, 
%                              measured by normalized decrease in variance 
%                              explained relative to the full model.
%             within.LL.raw -- (1 x numGroups) cell array; prominence of 
%                              dimGroups_within{i} evaluated on raw data, 
%                              measured by decrease in log-likelihood 
%                              relative to the full model.
%             within.VE.raw -- (1 x numGroups) cell array; prominence of
%                              dimGroups_within{i} evaluated on raw data, 
%                              measured by normalized decrease in variance
%                              explained relative to the full model.
%             (default: [])
%     metric -- string; relevant only if prom is specified. Specify which 
%               prominence metric to use for shading ('LL' or 'VE') 
%               (default: 'VE')
%
% Outputs:
%     
%    gp_params -- structure containing DLAG GP parameters, converted into
%                 units of time.
%                 DelayMatrix  -- (numGroups x xDim_across) array;
%                                 delays from across-group latents to 
%                                 observed variables
%                 tau_across -- (1 x xDim_across) array; across-group GP
%                               timescales
%                 tau_within -- (1 x numGroups) cell array; within-group
%                               GP timescales for each group. tau_within(i)
%                               is empty for groups with no within-group
%                               latents.
%              
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     10 Apr 2020 -- Initial full revision.
%     17 Apr 2020 -- Added 0-within-group dimension functionality.
%     28 Jun 2020 -- Added 0-across-group dimension functionality.
%     28 Sep 2020 -- Added option to shade latents according to prominence.

% Set optional arguments
plotAcross = true;
plotWithin = true;
units = '';
prom = [];
metric = 'VE';
assignopts(who,varargin);

xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
yDims = params.yDims;
numGroups = length(yDims);
colors = generateColors(); % Generate custom plotting colors
pointsize = 25; % Size of scatterplot points

% Convert GP params into units of time
gp_params = getGPparams_dlag(params, binWidth);

if xDim_across > 0
    % Take only the delays associated with the groups in rGroups
    delays = gp_params.DelayMatrix(rGroups(2),:) ...
           - gp_params.DelayMatrix(rGroups(1),:);

    % Figure out plot limits
    maxDelay = max(delays);
    minDelay = min(delays);
end

all_tau = [gp_params.tau_across gp_params.tau_within{:}];
maxTau = max(all_tau);

% Determine number of subplots, based on optional arguments
numPlot = sum([plotAcross*(xDim_across > 0) plotWithin.*(sum(xDim_within > 0))]);

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
if plotAcross && (xDim_across > 0)
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
    if isempty(prom)
        scatter(delays, gp_params.tau_across, pointsize, ...
                 'MarkerFaceColor', colors.grays{1}, ...
                 'MarkerEdgeColor', colors.grays{1});
    else
        % Shade latents according to prominence.
        c = prom.across.(metric).raw ./ max([prom.across.(metric).raw prom.within.(metric).raw{:}]);
        scatter(delays, gp_params.tau_across, pointsize, c, 'filled', 'MarkerEdgeColor', colors.grays{1});
        colormap(flipud(gray));
        colorbar;
        caxis([0 1]);
    end
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
            if isempty(prom)
                h = bar([1:xDim_within(groupIdx)],gp_params.tau_within{groupIdx},0.4);    
                set(h,'facecolor',colors.grays{1},'edgecolor',colors.grays{1});
            else
                for withinIdx = 1:xDim_within(groupIdx)
                    h = bar(withinIdx,gp_params.tau_within{groupIdx}(withinIdx),0.4);
                    set(h,'facecolor',colors.grays{1},'edgecolor',colors.grays{1});
                    set(h,'facealpha', ...
                        prom.within.(metric).raw{groupIdx}(withinIdx) / max([prom.across.(metric).raw prom.within.(metric).raw{:}])); 
                end
            end
            ylabel(sprintf('GP timescale%s', units));
            set(gca,'XTick',1:xDim_within(groupIdx));
            xlabel(sprintf('Within-group latents, group %d', groupIdx)); 
            hold off;
        end
    end
end