function plotDelaySignificance(sig, params, binWidth, rGroups)
%
% plotDelaySignificance(sig)
%
% Description: Plot the bootstrapped significance of the delays for each 
%              across-group latent.
%
% Arguments:
%
%     sig  -- structure containing the following fields:
%             raw    -- (1 x xDim_across) array; significance of each delay
%                       evaluated on raw data, measured by decrease in 
%                       log-likelihood relative to the unaltered model.
%             upper  -- (1 x xDim_across) array; upper bound of bootstrap 
%                       confidence interval
%             lower  -- (1 x xDim_across) array; lower bound of bootstrap 
%                       confidence interval
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
%     binWidth -- float; bin width or sample period (in e.g., ms)
%     rGroups  -- (1 x 2) array; Indexes of groups to get relative
%                 delays. rGroups(1) is the reference group.
%
% Outputs:
%
%     None. (But creates figures)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 May 2020 -- Initial full revision.

xDim_across = length(sig.raw);
colors = generateColors(); % Generate custom plotting colors

% Determine vertical axis limits for all subplots
upper = max(sig.upper);
lower = min([sig.lower 0]); % Include 0, if all lower bounds are above 0
range = upper - lower;
vlim = ([lower-0.1*range upper+0.1*range]);

figure;
hold on;

% Format subplot
xlim([0 xDim_across+1]);
set(gca, 'XTick', 1:xDim_across);

xtcklbl = cell(2,xDim_across);
% Get GP params for x tick labels
gp_params = getGPparams_dlag(params, binWidth);
% Take only the delays associated with the groups in rGroups
delays = gp_params.DelayMatrix(rGroups(2),:) ...
       - gp_params.DelayMatrix(rGroups(1),:);
for i = 1:xDim_across
    s1 = sprintf('%d: D=%.01f, ', i, delays(i));
    s2 = sprintf('=%.01f', gp_params.tau_across(i));
    xtcklbl{i} = [s1 '\tau' s2];
end
set(gca, 'XTickLabel', xtcklbl);

xlabel('Across-group latents');
ylim(vlim);
ylabel('Delay significance (\DeltaLL)');

% Plot data
line([0 xDim_across+1], [0 0], ...
     'Color', colors.grays{6}, ...
     'linestyle', '--', ...
     'linewidth', 1.5);
errorbar(1:xDim_across, ...
         sig.raw, ...
         sig.raw - sig.lower, ...
         sig.upper - sig.raw, ...
         'o', ...
         'Color', colors.grays{1}, ...
         'linestyle', 'none', ...
         'linewidth', 1.5, ...
         'markerfacecolor', colors.grays{1});
hold off;