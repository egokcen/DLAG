function plotSeqRaster(y, binWidth, varargin)
%
% plotSeqRaster(seq, binWidth, ...)
%
% Description: Plot all dimensions of a sequence in a raster-style plot.
%
% Arguments:
%
%     Required:
%
%     y        -- (yDim x T) array; sequence to be plotted
%     binWidth -- float; bin width or sample period, in units of time
%
%     Optional:
%
%     units    -- string; units of time of binWidth (for labels)
%                 (default: '')
%     
% Outputs:
%     None. (But creates a figure)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     11 Apr 2020 -- Initial full revision. 

units = '';
assignopts(who, varargin);

% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

[yDim, T] = size(y);
colors = generateColors(); % Get custom colors for plotting

% Set up x- and y-axis ticks
xtk = [1 T].*binWidth;
xtkl = [0 T-1].*binWidth;
ytk = 1:2:2*yDim;
ytkl = 1:yDim;

% Create raster-like data. Each row is zero-centered and normalized
% according to its own mean and absolute maximum.
Ynorm = (y - repmat(mean(y,2),1,T)) ./ max(abs(y),[],2) + ytk';

% Plot sequence rasters
figure;
h = subplot(1,1,1);
hold on;
plot((1:T).*binWidth, Ynorm(:,:), 'color', colors.grays{1}, 'linewidth', 1.0);
% Format axes
xlabel(sprintf('Time%s', units));
ylabel('Sequences');
set(h, 'xtick', xtk, 'xticklabel', xtkl, 'ytick', ytk, 'yticklabel', ytkl);