function zero_corr = zeroDelayCorr(ccmap, varargin)
%
% zero_corr = zeroDelayCorr(ccmap, ...)
%
% Description: Given a series of cross-correlation maps, compute the
%              zero-lag correlation at each time point.
%
% Arguments:
%
%     Required:
%
%     ccmap -- (1 x xDim_across) cell array; ccmat{i} -- (T x 2*T-1) array;
%              cross-correlation map between a pair of groups through
%              across-group latent i ("unwrapped" version of 'ccmat').
%
%     Optional:
%
%     binWidth -- float; bin width or sample period, in units of time
%                 (default: 1)
%     showPlot -- logical; set true to plot effective timescales. 
%                 (default: false)
%     units    -- string; if plotting, specify the units of time, e.g., 
%                 'ms' (default: [])
%
% Outputs:
%
%     zero_corr -- (1 x xDim) cell array; zero_corr{i} -- (1 x T) array; 
%                  zero-delay correlation at each time point for latent i.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2021 -- Initial full revision.

binWidth = 1;
showPlot = false;
units = [];
assignopts(who,varargin);

% Constants
xDim = length(ccmap);
T = size(ccmap{1},1); % Length of sequence (trial)

zero_corr = cell(1,xDim);
for j = 1:xDim
    curr_ccmap = flipud(ccmap{j});
    zero_corr{j} = curr_ccmap(:,T);
end

if showPlot
    zero_corr_all = [zero_corr{:}];
    ymax = 1.05.*max(zero_corr_all(:));
    figure;
    for j = 1:xDim
        subplot(1,xDim,j);
        hold on;
        plot((1:T).*binWidth, zero_corr{j}, 'k-');
        xlabel(sprintf('Time%s', [' (' units ')']));
        ylabel('Zero-lag correlation');
        ylim([0 ymax]);
    end
end