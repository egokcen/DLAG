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
%     28 Apr 2021 -- Added exception handling for xDim = 0 case.
%                    Improved plot formatting along vertical axis.
%     02 Feb 2022 -- Modified output units of delay_eff to match the units
%                    binWidth.

showPlot = false;
units = [];
assignopts(who,varargin);

% Constants
xDim = length(ccmap);
if xDim <= 0
   fprintf('computeEffectiveDelay: xDim = 0. Returning empty structure delay_eff\n');
   delay_eff = {};
   return;
end
T = size(ccmap{1},1); % Length of sequence (trial)

delay_eff = cell(1,xDim);
for j = 1:xDim
    [~, delay_eff{j}] = max(flipud(ccmap{j}),[],2);
    delay_eff{j} = (delay_eff{j}-T).*binWidth;
end

if showPlot
    % Determine vertical axis limits
    delay_all = [delay_eff{:}];
    ymax = 1.05.*max([abs(delay_all(:)); binWidth]);
    figure;
    for j = 1:xDim
        subplot(1,xDim,j);
        hold on;
        plot((1:T).*binWidth, delay_eff{j}, 'k-');
        plot((1:T).*binWidth,zeros(1,T),'color',[0.5 0.5 0.5],'linestyle','--');
        xlabel(sprintf('Time%s', [' (' units ')']));
        ylabel(sprintf('Effective Delay%s', [' (' units ')']));
        ylim([-ymax ymax]);
    end
end