function plotCrossCorrMap(ccmap, binWidth, varargin)
%
% plotCrossCorrMap(ccmap, binWidth, ...)
%
% Description: Plot a series of cross-correlation maps.
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
%     peak -- logical; set to true to superimpose a trace of the delay
%             at which the peak correlation occurs at each time point.
%             (default: false)
%     units -- string; specify the units of time, e.g., 'ms' (default: [])
%
% Outputs:
%
%     None. (Generates figures)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2021 -- Initial full revision.

units = [];
peak = false;
assignopts(who,varargin);

% Constants
xDim = length(ccmap);
T = size(ccmap{1},1); % Length of sequence (trial)
t0 = T;               % 0-delay index

figure;
for j = 1:xDim
    subplot(1,xDim,j);
    hold on;
    imagesc(flipud(ccmap{j}(:,:)));
    plot(t0.*ones(1,T),1:T,'color', [0.5 0.5 0.5], 'linestyle', '--');
    if peak
        % Compute delay at which peaks occur
        [~, peakdelay] = max(flipud(ccmap{j}),[],2);
        plot(peakdelay,1:T, 'k-');
    end
    colorbar;
    colormap parula;
    xticks([1 T 2*T-1]);
    xticklabels([-(T-1) 0 T-1].*binWidth);
    yticks([1 T]);
    yticklabels([0 T-1].*binWidth);
    xlabel(sprintf('Delay%s', [' (' units ')']));
    ylabel(sprintf('Time%s', [' (' units ')']));
end