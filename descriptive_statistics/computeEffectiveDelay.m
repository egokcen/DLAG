function delay_eff = computeEffectiveDelay(ccmap, binWidth, varargin)
%
% delay_eff = computeEffectiveDelay(ccmap, binWidth, ...)
%
% Description: Given a series of cross-correlation maps, compute the
%              effective delay at each time point.
%
% Arguments:
%
%     Required:
%
%     ccmap -- (1 x xDim_across) cell array; ccmat{i} -- (T x 2*T-1) array;
%              cross-correlation map between a pair of groups through
%              across-group latent i ("unwrapped" version of 'ccmat').
%     binWidth -- float; bin width or sample period, in units of time
%
%     Optional:
%
%     showPlot -- logical; set true to plot effective timescales. 
%                 (default: false)
%     units    -- string; if plotting, specify the units of time, e.g., 
%                 'ms' (default: [])
%
% Outputs:
%
%     delay_eff -- (1 x xDim) cell array; delay_eff{i} -- (1 x T) array; 
%                  effective delay at each time point for latent i.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2021 -- Initial full revision.

showPlot = false;
units = [];
assignopts(who,varargin);

% Constants
xDim = length(ccmap);
T = size(ccmap{1},1); % Length of sequence (trial)

delay_eff = cell(1,xDim);
for j = 1:xDim
    [~, delay_eff{j}] = max(flipud(ccmap{j}),[],2);
end

if showPlot
    figure;
    for j = 1:xDim
        subplot(1,xDim,j);
        hold on;
        plot((1:T).*binWidth, (delay_eff{j}-T).*binWidth, 'k-');
        plot((1:T).*binWidth,zeros(1,T),'color',[0.5 0.5 0.5],'linestyle','--');
        xlabel(sprintf('Time%s', [' (' units ')']));
        ylabel(sprintf('Effective Delay%s', [' (' units ')']));
        ymax = 1.05.*max([abs((delay_eff{j}-T).*binWidth); binWidth]);
        ylim([-ymax ymax]);
    end
end