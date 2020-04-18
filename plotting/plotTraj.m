function plotTraj(seq, xspec, varargin)
%
% plotTraj(seq, xspec, ...)
%
% Description: Plot neural trajectories in a two- or three-dimensional 
%              space.
%
% Arguments:
%
%     Required:
%
%     seq        -- data structure containing extracted trajectories
%     xspec      -- string; field name of trajectories in 'seq' to be 
%                   plotted (e.g., 'xorth' or 'xsm')
%
%     Optional:
%
%     dimsToPlot -- selects three dimensions in seq.(xspec) to plot 
%                   (default: 1:3)
%     nPlotMax   -- int; maximum number of trials to plot (default: 20)
%     plotSingle -- logical; if true, plot single-trial trajectories
%                   (default: true)
%     plotMean   -- logical; if true, plot mean across single-trial
%                   trajectories. Mean will be over the nPlotMax trials
%                   being plotted. (default: true)
%     plotZero  -- logical; if true, plot zero-firing rate point on each 
%                  plot (default: false)
%     Zero      -- (xDim x 1) array; The zero-firing rate point projected 
%                  into latent space (default: [])
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     11 Apr 2020 -- Initial full revision.

dimsToPlot = 1:3;
nPlotMax   = 20;
plotSingle = true;
plotMean   = true;
plotZero   = false;
Zero       = [];
assignopts(who, varargin);

if size(seq(1).(xspec), 1) < 2
    fprintf('ERROR: Trajectories have less than 2 dimensions.\n');
    return
end

colors = generateColors(); % Get custom colors for plotting

f = figure;
pos = get(gcf, 'position');
set(f, 'position', [pos(1) pos(2) 1.3*pos(3) 1.3*pos(4)]);
hold on;

Tmin = min([seq.T]);
% If trials are different lengths, compute mean over only the first Tmin 
% time points
xmean = zeros(length(dimsToPlot),Tmin); 
N = min(length(seq), nPlotMax);
if plotZero
    % Plot the zero-firing rate point in latent space
    Zero = Zero(dimsToPlot,:);
    if length(dimsToPlot) <= 2
        plot(Zero(1), Zero(2), 'd', 'Color', colors.grays{1}, ...
             'MarkerSize', 10.0, 'MarkerFaceColor', colors.grays{1});
    else
        plot3(Zero(1), Zero(2), Zero(3), 'd', 'Color', colors.grays{1}, ...
             'MarkerSize', 10.0, 'MarkerFaceColor', colors.grays{1});
    end
end
for n = 1:N
    dat = seq(n).(xspec)(dimsToPlot,:);
    if plotSingle
        % Plot single-trial trajectories
        if length(dimsToPlot) <=2
            plot(dat(1,:), dat(2,:), '.-', ...
                 'linewidth', 0.5, ...
                 'color', colors.blues{7});
        else
            plot3(dat(1,:), dat(2,:), dat(3,:), '.-', ...
                  'linewidth', 0.5, ...
                  'color', colors.blues{7});
        end
    end
    for i = 1:length(dimsToPlot)
        xmean(i,:) = xmean(i,:) + (1.0/N) .* dat(i,1:Tmin); 
    end
end
if plotMean
    % Plot the mean trajectory
    if length(dimsToPlot) <= 2
        plot(xmean(1,:), xmean(2,:), '.-', ...
             'MarkerSize', 20.0, ...
             'linewidth', 2.0, ...
             'color', colors.blues{4}); 
    else
        plot3(xmean(1,:), xmean(2,:), xmean(3,:), '.-', ...
              'MarkerSize', 20.0, ...
              'linewidth', 2.0, ...
              'color', colors.blues{4});
    end
end

% Additional formatting of axis labels.
axis equal;
grid on;
str1 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(1));
str2 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(2));
xlabel(str1, 'interpreter', 'latex', 'fontsize', 24);
ylabel(str2, 'interpreter', 'latex', 'fontsize', 24);
if length(dimsToPlot) > 2
    str3 = sprintf('$${\\mathbf x}_{%d,:}$$', dimsToPlot(3));
    zlabel(str3, 'interpreter', 'latex', 'fontsize', 24);
end