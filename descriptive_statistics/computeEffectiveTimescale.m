function tau_eff = computeEffectiveTimescale(ccmap, binWidth, varargin)
%
% tau_eff = computeEffectiveTimescale(ccmap, binWidth, ...)
%
% Description: Given a series of cross-correlation maps, compute the
%              effective timescale at each time point.
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
%     tau_eff -- (1 x xDim) cell array; tau_eff{i} -- (1 x T) array; 
%                effective timescale at each time point for latent i.
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

tau_eff = cell(1,xDim);
for j = 1:xDim
    tau_eff{j} = nan(1,T);
    for t = 1:T
        % Get a slice of the cross-correlation map
        ccmap_ud = flipud(ccmap{j}); % The ccg is upside down wrt time.
        ccmap_slice = ccmap_ud(t,:);
        % Find the max correlation
        [maxCorr, maxIdx] = max(ccmap_slice);
        % Find all time points less than 1/sqrt(e) of the max correlation
        tIdxs = find(ccmap_slice <= maxCorr/sqrt(exp(1)));
        % Find time differences between these points and the max point
        tDiffs = abs(maxIdx - tIdxs);
        % Take the effective timescale to be the difference to the time
        % point closest to the max correlation time point.
        if isempty(tDiffs)
            tau_eff{j}(t) = nan;
        else
            tau_eff{j}(t) = min(tDiffs).*binWidth;
        end
    end
end

if showPlot
    figure;
    for j = 1:xDim
        subplot(1,xDim,j);
        hold on;
        plot((1:T).*binWidth, tau_eff{j}(1:T),'k-');
        xlabel(sprintf('Time%s', [' (' units ')']));
        ylabel(sprintf('Effective GP Timescale%s', [' (' units ')']));
        ylim([0 max(tau_eff{j})]);
    end
end