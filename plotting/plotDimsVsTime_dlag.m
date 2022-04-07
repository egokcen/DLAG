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
%                  NOTE: Not relevant if trialGroups is specified.
%     plotSingle -- logical; if true, plot single-trial trajectories
%                   (default: true)
%     plotMean   -- logical; if true, plot mean across single-trial
%                   trajectories. Mean will be over the nPlotMax trials
%                   being plotted. (default: true)
%     units     -- string; units of time of binWidth (for labels)
%                  (default: '')
%     trialGroups -- (1 x numTrialGroups) cell array; Each element contains
%                    a list of trials to be grouped together for
%                    color-coding and for computing a mean time course.
%                    (default: {})
%     trialColors -- (1 x numTrialGroups) cell array; If trialGroups is
%                    specified, then trialColors must be specified, 
%                    where each element is the color for a given trial
%                    group. (default: {})
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     19 Mar 2021 -- Initial full revision.
%     10 Apr 2021 -- Added redTrials optional argument.
%     11 Apr 2021 -- Fixed issue with plot gap when xDim_across == 0.
%     12 Apr 2021 -- Modified position of figure (issues on some OS's)
%     13 Sep 2021 -- Updated labeling conventions.
%     31 Mar 2022 -- Added option to group and color-code trials.

nPlotMax   = 20;
plotSingle = true;
plotMean   = true;
units      = '';
trialGroups = {};
trialColors = {};
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
N = 0;
allTrials = [];
if isempty(trialGroups)
    N = min(length(seq), nPlotMax); 
    allTrials = 1:N;
else
    for trialGroupIdx = 1:length(trialGroups)
        N = N + length(trialGroups{trialGroupIdx}); 
        allTrials = [allTrials trialGroups{trialGroupIdx}];
    end
end
  
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
    'OuterPosition', [0.05 0.05 min([1 0.2*nCols]) min([1 0.25*nRows])]);
Tmax    = max([seq.T]);  
Tmin    = min([seq.T]);
xtkStep = ceil(Tmax/25)*5;
xtk     = 1:xtkStep:Tmax;
xtkl    = 0:(xtkStep*binWidth):(Tmax-1)*binWidth;

% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

fontsize = 12; % Size of text label fonts

for groupIdx = 1:numGroups

    % Determine vertical scale for current group
    Xall = [groupSeq{groupIdx}(allTrials).(xspec)];
    xMax = ceil(10*max(abs(Xall(:)))) / 10;
    xMid = xMax / 2; %ceil(10*(xMax/2)) / 10;
    ytk     = [-xMax -xMid 0 xMid xMax];

    % Convert delay labels to units of time
    gp_params = getGPparams_dlag(groupParams{groupIdx}, binWidth);
    DelayMatrix = gp_params.DelayMatrix;

    % Partition latents into across- and within-group components
    [seqAcross, seqWithin] = partitionLatents_meanOnly(groupSeq{groupIdx}, ...
        xDim_across, xDim_within(groupIdx), 'xspec', xspec);

    % Initialize plotIdx to first column of the current row
    plotIdx = (groupIdx-1)*nCols;
    
    % Across-group latents
    for k = 1:xDim_across
        plotIdx = plotIdx + 1; % Increment to the next column
        h = subplot(nRows, nCols, plotIdx);
        hold on;
        if isempty(trialGroups)
            % Initialize the mean trajectory. Only average over time points up to
            % Tmin, if trial lengths are different.
            xmean = zeros(1,Tmin); 
            for n = allTrials
                dat = seqAcross(n).(xspec);
                if plotSingle 
                    % Plot single-trial trajectories
                    T = seqAcross(n).T;
                    col = colors.grays{5};
                    plot(1:T, dat(k,:), ...
                         'linewidth', 0.05, ...
                         'color', col);
                end
                xmean = xmean + (1.0/N) .* dat(k,1:Tmin);
            end
            % Plot the mean trajectory
            if plotMean
                plot(1:Tmin, xmean, ...
                     'linewidth', 2.0, ... 
                     'color', colors.grays{1});
            end     
        else
            for trialGroupIdx = 1:length(trialGroups)
                % Initialize the mean trajectory. Only average over time 
                % points up to Tmin, if trial lengths are different.
                xmean = zeros(1,Tmin); 
                for n = trialGroups{trialGroupIdx}
                    dat = seqAcross(n).(xspec);
                    if plotSingle 
                        % Plot single-trial trajectories
                        T = seqAcross(n).T;
                        plot(1:T, dat(k,:), ...
                             'linewidth', 0.05, ...
                             'color', trialColors{trialGroupIdx});
                    end
                    xmean = xmean + (1.0/length(trialGroups{trialGroupIdx})) ...
                            .* dat(k,1:Tmin);
                end
                % Plot the mean trajectory
                if plotMean
                    plot(1:Tmin, xmean, ...
                         'linewidth', 2.0, ... 
                         'color', trialColors{trialGroupIdx});
                end
                
            end
        end
        % Additional formatting of titles and axis labels.
        axis([1 Tmax 1.1*min(ytk) 1.1*max(ytk)]);
        set(h, 'xtick', xtk, 'xticklabel', xtkl);
        set(h, 'ytick', ytk, 'yticklabel', ytk);
        str = sprintf('$${\\mathbf x}^{a}_{%d,%d,:}$$',groupIdx,k);
        str = sprintf('%s, $$D_{%d,%d} = %3.1f$$', str, groupIdx, k, DelayMatrix(k));
        str = sprintf('%s%s', str, units);
        title(str, 'interpreter', 'latex', 'fontsize', fontsize);
        xlabel(sprintf('Time%s', units));
    end
    
    % Leave a one-panel gap between across- and within-group latents
    if xDim_across > 0
        plotIdx = plotIdx+1;
    end
    
    % Within-group latents
    for k = 1:xDim_within(groupIdx)
        plotIdx = plotIdx + 1;
        h = subplot(nRows, nCols, plotIdx);
        hold on;
        if isempty(trialGroups)
            % Initialize the mean trajectory. Only average over time points up to
            % Tmin, if trial lengths are different.
            xmean = zeros(1,Tmin); 
            for n = allTrials
                dat = seqWithin{1}(n).(xspec);
                if plotSingle 
                    % Plot single-trial trajectories
                    T = seqWithin{1}(n).T;
                    col = colors.grays{5};
                    plot(1:T, dat(k,:), ...
                         'linewidth', 0.05, ...
                         'color', col);
                end
                xmean = xmean + (1.0/N) .* dat(k,1:Tmin);
            end
            % Plot the mean trajectory
            if plotMean
                plot(1:Tmin, xmean, ...
                     'linewidth', 2.0, ... 
                     'color', colors.grays{1});
            end
            
        else
            for trialGroupIdx = 1:length(trialGroups)
                % Initialize the mean trajectory. Only average over time 
                % points up to Tmin, if trial lengths are different.
                xmean = zeros(1,Tmin); 
                for n = trialGroups{trialGroupIdx}
                    dat = seqWithin{1}(n).(xspec);
                    if plotSingle 
                        % Plot single-trial trajectories
                        T = seqWithin{1}(n).T;
                        plot(1:T, dat(k,:), ...
                             'linewidth', 0.05, ...
                             'color', trialColors{trialGroupIdx});
                    end
                    xmean = xmean + (1.0/length(trialGroups{trialGroupIdx})) ...
                            .* dat(k,1:Tmin);
                end
                % Plot the mean trajectory
                if plotMean
                    plot(1:Tmin, xmean, ...
                         'linewidth', 2.0, ... 
                         'color', trialColors{trialGroupIdx});
                end 
            end
        end
        
        % Additional formatting of titles and axis labels.
        axis([1 Tmax 1.1*min(ytk) 1.1*max(ytk)]);
        set(h, 'xtick', xtk, 'xticklabel', xtkl);
        set(h, 'ytick', ytk, 'yticklabel', ytk);
        str = sprintf('$${\\mathbf x}^{w}_{%d,%d,:}$$',groupIdx,k);
        title(str, 'interpreter', 'latex', 'fontsize', fontsize);
        xlabel(sprintf('Time%s', units));
    end

end