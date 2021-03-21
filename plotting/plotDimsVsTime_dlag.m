function plotDimsVsTime_dlag(seq, xspec, params, binWidth, varargin)
%
% plotDimsVsTime_dlag(seq, xspec, binWidth, varargin)
%
% Description: Plot each DLAG dimension versus time in a separate panel, 
%              along with the mean trajectory across trials.
%              Group and scale panels according to which observation group
%              trajectories belong to.
%
% Arguments:
%
%     Required:
%
%     seq       -- data structure containing extracted trajectories
%     xspec     -- string; field name of trajectories in 'seq' to be 
%                  plotted (e.g., 'xorth' or 'xsm')
%     binWidth  -- bin width used when fitting model
%     params    -- Structure containing DLAG model parameters.
%                    DelayMatrix  -- (numGroups x xDim_across) array;
%                                    delays from across-group latents to 
%                                    observed variables. NOTE: Delays are
%                                    reported as (real-valued) number of
%                                    time-steps.
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number of
%                                    within-group latents in each group
%
%     Optional:
%
%     nPlotMax  -- int; maximum number of trials to plot (default: 20)
%     plotSingle -- logical; if true, plot single-trial trajectories
%                   (default: true)
%     plotMean   -- logical; if true, plot mean across single-trial
%                   trajectories. Mean will be over the nPlotMax trials
%                   being plotted. (default: true)
%     units     -- string; units of time of binWidth (for labels)
%                  (default: '')
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     19 Mar 2021 -- Initial full revision.

nPlotMax   = 20;
plotSingle = true;
plotMean   = true;
units      = '';
assignopts(who, varargin);

colors = generateColors(); % Get custom colors for plotting
numGroups = length(params.yDims);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
xDim_total = xDim_across + xDim_within;

% Check if there are any trajectories
Xall = [seq.(xspec)];

% Set up figure and axes
if isempty(Xall)
    fprintf('plotDimsVsTime_dlag: No trajectories to plot.\n');
    return;
end

% Number of trajectories to be plotted
N = min(length(seq), nPlotMax); 
  
% Group sequences and parameters
groupSeq = partitionObs(seq,xDim_total,'datafield',xspec);
groupParams = partitionParams_dlag(params);

% Global figure properties
f = figure;

% Number of rows and columns
nRows = numGroups;
nCols = max(xDim_total);
if max(xDim_within) > 0 && xDim_across > 0
    % Put a gap between across- and within-group, assuming both types exist
    nCols = nCols + 1; 
end

% Size and axis scales
set(f, 'Units', 'Normalized', ...
    'OuterPosition', [1 1 min([1 0.2*nCols]) min([1 0.25*nRows])]);
Tmax    = max([seq.T]);  
Tmin    = min([seq.T]);
xtkStep = ceil(Tmax/25)*5;
xtk     = 1:xtkStep:Tmax;
xtkl    = 0:(xtkStep*binWidth):(Tmax-1)*binWidth;

% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

for groupIdx = 1:numGroups

    % Determine vertical scale for current group
    Xall = [groupSeq{groupIdx}(1:N).(xspec)];
    xMax = ceil(10*max(abs(Xall(:)))) / 10;
    xMid = ceil(10*(xMax/2)) / 10;
    ytk     = [-xMax -xMid 0 xMid xMax];

    % Convert delay labels to units of time
    gp_params = getGPparams_dlag(groupParams{groupIdx}, binWidth);
    DelayMatrix = gp_params.DelayMatrix;

    % Partition latents into across- and within-group components
    [seqAcross, seqWithin] = partitionLatents_meanOnly(groupSeq{groupIdx}(1:N), ...
        xDim_across, xDim_within(groupIdx), 'xspec', xspec);

    % Initialize plotIdx to first column of the current row
    plotIdx = (groupIdx-1)*nCols;
    
    % Across-group latents
    for k = 1:xDim_across
        plotIdx = plotIdx + 1; % Increment to the next column
        h = subplot(nRows, nCols, plotIdx);
        hold on;
        % Initialize the mean trajectory. Only average over time points up to
        % Tmin, if trial lengths are different.
        xmean = zeros(1,Tmin); 
        for n = 1:N
            dat = seqAcross(n).(xspec);
            if plotSingle 
                % Plot single-trial trajectories
                T = seqAcross(n).T;
                plot(1:T, dat(k,:), ...
                     'linewidth', 0.05, ...
                     'color', colors.grays{5});
            end
            xmean = xmean + (1.0/N) .* dat(k,1:Tmin);
        end
        % Plot the mean trajectory
        if plotMean
            plot(1:Tmin, xmean, ...
                 'linewidth', 2.0, ... 
                 'color', colors.grays{1});
        end

        % Additional formatting of titles and axis labels.
        axis([1 Tmax 1.1*min(ytk) 1.1*max(ytk)]);
        set(h, 'xtick', xtk, 'xticklabel', xtkl);
        set(h, 'ytick', ytk, 'yticklabel', ytk);
        str = sprintf('$${\\mathbf x}^{(%d,%d)}_{a,:}$$',groupIdx,k);
        str = sprintf('%s, $$D_{%d%d} = %3.1f$$', str, groupIdx, k, DelayMatrix(k));
        str = sprintf('%s%s', str, units);
        title(str, 'interpreter', 'latex', 'fontsize', 16);
        xlabel(sprintf('Time%s', units));
    end
    
    % Leave a one-panel gap between across- and within-group latents
    plotIdx = plotIdx+1;

    % Within-group latents
    for k = 1:xDim_within(groupIdx)
        plotIdx = plotIdx + 1;
        h = subplot(nRows, nCols, plotIdx);
        hold on;
        % Initialize the mean trajectory. Only average over time points up to
        % Tmin, if trial lengths are different.
        xmean = zeros(1,Tmin); 
        for n = 1:N
            dat = seqWithin{1}(n).(xspec);
            if plotSingle 
                % Plot single-trial trajectories
                T = seqWithin{1}(n).T;
                plot(1:T, dat(k,:), ...
                     'linewidth', 0.05, ...
                     'color', colors.grays{5});
            end
            xmean = xmean + (1.0/N) .* dat(k,1:Tmin);
        end
        % Plot the mean trajectory
        if plotMean
            plot(1:Tmin, xmean, ...
                 'linewidth', 2.0, ... 
                 'color', colors.grays{1});
        end

        % Additional formatting of titles and axis labels.
        axis([1 Tmax 1.1*min(ytk) 1.1*max(ytk)]);
        set(h, 'xtick', xtk, 'xticklabel', xtkl);
        set(h, 'ytick', ytk, 'yticklabel', ytk);
        str = sprintf('$${\\mathbf x}^{(%d,%d)}_{w,:}$$',groupIdx,k);
        title(str, 'interpreter', 'latex', 'fontsize', 16);
        xlabel(sprintf('Time%s', units));
    end

end