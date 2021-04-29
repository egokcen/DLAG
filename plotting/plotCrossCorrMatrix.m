function plotCrossCorrMatrix(ccmat, binWidth, varargin)
%
% plotCrossCorrMatrix(ccmat, binWidth, ...)
%
% Description: Plot a series of cross-correlation matrices.
%
% Arguments:
%
%     Required:
%
%     ccmat -- (1 x xDim_across) cell array; ccmat{i} -- (T x T) array;
%              cross-correlation matrix between a pair of groups through
%              across-group latent i.
%     binWidth -- float; bin width or sample period, in units of time
%
%     Optional:
%
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
%     28 Apr 2021 -- Added exception handling for xDim = 0 case.

units = [];
assignopts(who,varargin);

% Constants
xDim = length(ccmat);
if xDim <= 0
   fprintf('plotCrossCorrMatrix: xDim = 0. No cross-correlation matrices to plot.\n');
   return;
end
T = size(ccmat{1},1); % Length of sequence (trial)

figure;
for j = 1:xDim
    subplot(1,xDim,j);    
    hold on;
    imagesc(flipud(ccmat{j}));
    plot(1:T,fliplr(1:T),'linestyle', '--', 'color', 'k', 'linewidth', 1.0);
    colorbar;
    colormap parula;
    xticks([1 T]);
    xticklabels([0 T-1].*binWidth);
    yticks([1 T]);
    yticklabels([T-1 0].*binWidth);
    xlabel(sprintf('Time%s', [' (' units ')']));
    ylabel(sprintf('Time%s', [' (' units ')']));
end