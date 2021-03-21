function plotDimsVsTime(seq, xspec, binWidth, varargin)
%
% plotDimsVsTime(seq, xspec, binWidth, ...)
%
% Description: Plot generic (not DLAG-specific) state dimensions versus 
%              time in a separate panel, along with the mean trajectory 
%              across trials.
%
% Arguments:
%
%     Required:
%
%     seq       -- data structure containing extracted trajectories
%     xspec     -- string; field name of trajectories in 'seq' to be 
%                  plotted (e.g., 'xorth' or 'xsm')
%     binWidth  -- bin width used when fitting model
%
%     Optional:
%
%     nPlotMax  -- int; maximum number of trials to plot (default: 20)
%     nCols     -- int; number of subplot columns (default: 4)
%     plotSingle -- logical; if true, plot single-trial trajectories
%                   (default: true)
%     plotMean   -- logical; if true, plot mean across single-trial
%                   trajectories. Mean will be over the nPlotMax trials
%                   being plotted. (default: true)
%     plotZero  -- logical; if true, plot zero-firing rate point on each 
%                  plot (default: false)
%     Zero      -- (xDim x 1) array; The zero-firing rate point projected 
%                  into latent space (default: [])
%     units     -- string; units of time of binWidth (for labels)
%                  (default: '')
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     26 Mar 2020 -- Initial full revision.
%     17 Apr 2020 -- Added 0-within-group dimension functionality
%     28 Jun 2020 -- Added 0-across-group dimension functionality
%     21 Sep 2020 -- Fixed bug where across-group latents for first
%                    group only were plotted for all groups.
%     18 Mar 2021 -- Removed 'dlagFormat'. Updated some plotting 
%                    conventions.

nPlotMax   = 20;
nCols      = 4;
plotSingle = true;
plotMean   = true;
plotZero   = false;
Zero       = [];
units      = '';
assignopts(who, varargin);

colors = generateColors(); % Get custom colors for plotting

% Check if there are any trajectories
Xall = [seq.(xspec)];

% Set up figure and axes
if isempty(Xall)
    fprintf('plotDimsVsTime: No trajectories to plot.\n');
    return;
end

% Determine plot limits based only on trajectories that actually appear in
% the plot.
N = min(length(seq), nPlotMax); % Number of trajectories to be plotted
Xall = [seq(1:N).(xspec)];

% Determine number of subplot rows
nRows   = ceil(size(Xall, 1) / nCols);
    
f = figure;
% Set final figure size
set(f, 'Units', 'Normalized', ...
    'OuterPosition', [1 1 min([1 0.2*nCols]) min([1 0.25*nRows])]);
xMax = ceil(10*max(abs(Xall(:)))) / 10;
xMid = ceil(10*(xMax/2)) / 10;

Tmax    = max([seq.T]);  
Tmin    = min([seq.T]);
xtkStep = ceil(Tmax/25)*5;
xtk     = 1:xtkStep:Tmax;
xtkl    = 0:(xtkStep*binWidth):(Tmax-1)*binWidth;
ytk     = [-xMax -xMid 0 xMid xMax];

% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

for k = 1:size(Xall,1)
    subplot(nRows, nCols, k);
    hold on;
    if plotZero
        % Plot the zero-firing rate point in latent space
        plot([1 Tmin], [Zero(k) Zero(k)], ...
             'linestyle', '--', ...
             'color', colors.grays{1}, ...
             'linewidth', 1.5);
    end
    % Initialize the mean trajectory. Only average over time points up to
    % Tmin, if trial lengths are different.
    xmean = zeros(1,Tmin); 
    for n = 1:N
        dat = seq(n).(xspec);
        if plotSingle 
            % Plot single-trial trajectories
            T = seq(n).T;
            plot(1:T, dat(k,:), ...
                 'linewidth', 0.05, ...
                 'color', colors.grays{5});
        end
        xmean = xmean + (1.0/N) .* dat(k,1:Tmin);
    end
    if plotMean
        % Plot the mean trajectory
        plot(1:Tmin, xmean, ...
             'linewidth', 2.0, ... 
             'color', colors.grays{1});
    end
end

% Additional formatting of titles and axis labels.
for k = 1:size(Xall,1)
    h = subplot(nRows, nCols, k);
    axis([1 Tmax 1.1*min(ytk) 1.1*max(ytk)]);
    set(h, 'xtick', xtk, 'xticklabel', xtkl);
    set(h, 'ytick', ytk, 'yticklabel', ytk);
    str = sprintf('$${\\mathbf x}_{%d,:}$$',k);
    title(str, 'interpreter', 'latex', 'fontsize', 16);
    xlabel(sprintf('Time%s',units));
end
    