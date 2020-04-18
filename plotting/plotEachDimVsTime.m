function plotEachDimVsTime(seq, xspec, binWidth, varargin)
%
% plotEachDimVsTime(seq, xspec, binWidth, ...)
%
% Description: Plot each state dimension versus time in a separate panel, 
%              along with the mean trajectory across trials.
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
%     dlagFormat -- logical; if false, latents will plotted in the order as
%                   they appear in seq.(xspec), without any special 
%                   treatment. If true, figure will have the following
%                   format:
%                       row 1: across-group latents, group 1
%                       row 2: across-group latents, group 2
%                         .
%                         .
%                         .
%                       row numGroups+1: within-group latents, group 1
%                       row numGroups+2: within-group latents, group 2
%                         .
%                         .
%                         .
%                  (default: false)
%    params -- Structure containing DLAG model parameters. NOTE: If 
%              dlagFormat is set to true, then this argument needs to be
%              specified. Relevant fields:
% 
%                    DelayMatrix  -- (numGroups x xDim_across) array;
%                                    delays from across-group latents to 
%                                    observed variables. NOTE: Delays are
%                                    reported as (real-valued) number of
%                                    time-steps.
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number of
%                                    within-group latents in each group
%              (default: [])
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     26 Mar 2020 -- Initial full revision.
%     17 Apr 2020 -- Added 0-within-group dimension functionality


nPlotMax   = 20;
nCols      = 4;
plotSingle = true;
plotMean   = true;
plotZero   = false;
Zero       = [];
units      = '';
dlagFormat = false;
params     = [];
assignopts(who, varargin);

colors = generateColors(); % Get custom colors for plotting

Xall = [seq.(xspec)];

% Set up figure and axes
if isempty(Xall)
    fprintf('plotEachDimVsTime: No trajectories to plot.\n');
    return;
end
    
f = figure;
pos = get(gcf, 'position');
set(f, 'position', [pos(1) pos(2) 2*pos(3) pos(4)]);

xMax = ceil(10 * max(abs(Xall(:)))) / 10; % round max value to next highest 1e-1

Tmax    = max([seq.T]);  
Tmin    = min([seq.T]);
xtkStep = ceil(Tmax/25)*5;
xtk     = 1:xtkStep:Tmax;
xtkl    = 0:(xtkStep*binWidth):(Tmax-1)*binWidth;
ytk     = [-xMax 0 xMax];

N = min(length(seq), nPlotMax); % Number of trajectories to be plotted

% Format units for axis labels
if ~isempty(units)
    units = sprintf(' (%s)', units); 
end

if ~dlagFormat % Plot all trajectories without any special treatment
    % Determine number of subplot rows
    nRows   = ceil(size(Xall, 1) / nCols);
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
                     'color', colors.blues{7});
            end
            xmean = xmean + (1.0/N) .* dat(k,1:Tmin);
        end
        if plotMean
            % Plot the mean trajectory
            plot(1:Tmin, xmean, ...
                 'linewidth', 2.0, ... 
                 'color', colors.blues{4});
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
    
else % Format figure to more clearly distinguish across- and within-group latents
    
    assert(~isempty(params));
    xDim_across = params.xDim_across;
    xDim_within = params.xDim_within;
    % Convert delay labels to units of time
    gp_params = getGPparams_dlag(params, binWidth);
    DelayMatrix = gp_params.DelayMatrix;
    
    % First partition latents into across- and within-group components
    [seqAcross, seqWithin] = partitionLatents_meanOnly(seq, ...
        xDim_across, xDim_within, 'xspec', xspec);
   
    % Determine number of subplot rows and columns
    numGroups = length(xDim_within);
    nRows = 2 * numGroups;
    nCols = max([xDim_across xDim_within]);
    
    for groupIdx = 1:numGroups
        % Across-group latents
        for k = 1:xDim_across
            plotIdx = (groupIdx-1)*nCols + k;
            h = subplot(nRows, nCols, plotIdx);
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
                dat = seqAcross(n).(xspec);
                if plotSingle 
                    % Plot single-trial trajectories
                    T = seqAcross(n).T;
                    plot(1:T, dat(k,:), ...
                         'linewidth', 0.05, ...
                         'color', colors.blues{7});
                end
                xmean = xmean + (1.0/N) .* dat(k,1:Tmin);
            end
            % Plot the mean trajectory
            if plotMean
                plot(1:Tmin, xmean, ...
                     'linewidth', 2.0, ... 
                     'color', colors.blues{4});
            end
            
            % Additional formatting of titles and axis labels.
            axis([1 Tmax 1.1*min(ytk) 1.1*max(ytk)]);
            set(h, 'xtick', xtk, 'xticklabel', xtkl);
            set(h, 'ytick', ytk, 'yticklabel', ytk);
            str = sprintf('$${\\mathbf x}^{(%d,%d)}_{a,:}$$',groupIdx,k);
            str = sprintf('%s, $$D_{%d%d} = %3.1f$$', str, groupIdx, k, DelayMatrix(groupIdx,k));
            str = sprintf('%s%s', str, units);
            title(str, 'interpreter', 'latex', 'fontsize', 16);
            xlabel(sprintf('Time%s', units));
        end
        
        % If this group has no within-group latents, then leave this row
        % empty.
        if xDim_within(groupIdx) > 0
            % Within-group latents
            for k = 1:xDim_within(groupIdx)
                plotIdx = numGroups*xDim_across + (groupIdx-1)*nCols + k;
                h = subplot(nRows, nCols, plotIdx);
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
                    dat = seqWithin{1,groupIdx}(n).(xspec);
                    if plotSingle 
                        % Plot single-trial trajectories
                        T = seqWithin{1,groupIdx}(n).T;
                        plot(1:T, dat(k,:), ...
                             'linewidth', 0.05, ...
                             'color', colors.blues{7});
                    end
                    xmean = xmean + (1.0/N) .* dat(k,1:Tmin);
                end
                % Plot the mean trajectory
                if plotMean
                    plot(1:Tmin, xmean, ...
                         'linewidth', 2.0, ... 
                         'color', colors.blues{4});
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

    end
end