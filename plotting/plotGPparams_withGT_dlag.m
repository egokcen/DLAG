function plotGPparams_withGT_dlag(estParams,trueParams,binWidth,rGroups,varargin)
%
% plotGPparams_withGT_dlag(estParams,trueParams,binWidth,rGroups,...)
%
% Description: Plot DLAG Delay Matrix and across-group latent timescales
%              for both estimated and ground-truth parameters, on the same
%              plot. Within-group GP timescales are not plotted. Since
%              latent variables are not ordered, in general, plotting
%              within-group parameters on the same plot may only be
%              confusing.
%
% Arguments: 
%
%     Required:
%
%     estParams   -- Structure containing estimated DLAG model parameters. 
%                    Contains the fields
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
%                    xDim_within  -- (1 x numGroups) array; number 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
%     trueParams -- Structure containing ground truth DLAG model parameters.
%                   trueParams has the same format as estParams.
%
%     binWidth   -- float; bin width or sample period (in e.g., ms)
%     rGroups    -- (1 x 2) array; Indexes of groups to get relative
%                   delays. rGroups(1) is the reference group.
%
%     Optional:
%
%     units      -- string; units of time of binWidth (for labels)
%                   (default: '')
%
% Outputs:
%     None. (But creates figures)
%              
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     10 Apr 2020 -- Initial full revision.

% Set optional arguments
units = '';
assignopts(who,varargin);

colors = generateColors(); % Generate custom plotting colors

% Convert GP params into units of time
estParams = getGPparams_dlag(estParams, binWidth);
trueParams = getGPparams_dlag(trueParams, binWidth);

% Take only the delays associated with the groups in rGroups
delays_est = estParams.DelayMatrix(rGroups(2),:) ...
           - estParams.DelayMatrix(rGroups(1),:);
       
delays_true = trueParams.DelayMatrix(rGroups(2),:) ...
            - trueParams.DelayMatrix(rGroups(1),:);

% Figure out plot limits
maxDelay = max([delays_est delays_true]);
minDelay = min([delays_est delays_true]);

maxTau = max([estParams.tau_across trueParams.tau_across]);
minTau = min([estParams.tau_across trueParams.tau_across]);

figure;
hold on;

% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

xlabel(sprintf('Delay from area %d to area %d%s',rGroups(1),rGroups(2),units));
ylabel(sprintf('Across-group GP timescale%s', units));
xlim([minDelay-10,maxDelay+10]);
ylim([minTau-10,maxTau+10]);
line([0 0], [minTau-10 maxTau+10], ...
     'Color', colors.grays{6}, ...
     'linestyle', '--', ...
     'linewidth', 1.5);
% Ground truth parameters
h1 = scatter(delays_true, trueParams.tau_across, ...
             'MarkerFaceColor', colors.grays{1}, ...
             'MarkerEdgeColor', colors.grays{1});
% Estimated parameters
h2 = scatter(delays_est, estParams.tau_across, ...
             'MarkerFaceColor', colors.reds{3}, ...
             'MarkerEdgeColor', colors.reds{3});
legend([h1 h2], 'Ground-truth', 'Estimated');
hold off;
