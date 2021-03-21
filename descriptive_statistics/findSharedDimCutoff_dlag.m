function d_shared = findSharedDimCutoff_dlag(params, cutoffPC, varargin)
%
% d_shared = findSharedDimCutoff_dlag(params, cutoffPC, ...)
%
% Description: Find the number of overall dimensions (within or across) 
%              that explain a certain percentage of shared variance, 
%              for each group individually.
%              Also, find the numer of across-group dimensions that explain
%              a certain percentage of the shared across-group covariance.
%
% Arguments:
% 
%     Required:
%
%     params   -- Structure containing DLAG model parameters.
%                 Contains the (relevant) fields
% 
%                    C            -- (yDim x (numGroups*xDim)) array;
%                                    mapping between low- and high-d spaces
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number of 
%                                    within-group latents in each group
%     cutoffPC -- float; cutoff proportion of shared variance explained 
%                 (between 0 and 1)
% 
%     Optional:
%
%     plotSpec -- logical; if true, plot the spectra and cumulative shared
%                 variance explained. (default: false)
%     rGroups  -- (1 x 2) int array; Indexes of the two groups whose cross-
%                 covariance we wish to consider. (default: [1 2])
%
% Outputs:
%
%     d_shared -- structure with the following fields:
%                 total  -- (1 x numGroups) array; number of dimensions
%                            required to explain cutoffPC of the
%                            shared variance for each group individually.
%                 across -- int; number of dimensions required to explain
%                           cutoffPC of the across-group-covariance among
%                           the two groups specified in rGroups.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     11 Apr 2020 -- Initial full revision.
%     17 Apr 2020 -- Added 0-within-group dimension functionality
%     24 May 2020 -- Updated argument documentation.
%     17 Jun 2020 -- Major revision to compute cutoffs for overall
%                    dimensions, rather than within- or across-group only.
%                    Within- and across-group only cutoffs are poorly
%                    defined, since within- and across-group dimensions
%                    are, in general, correlated.
%     28 Jun 2020 -- Fixed issue with reporting cutoff dimensionality.
%     16 Aug 2020 -- Added d_shared.across functionality.
%     19 Mar 2021 -- Removed d_shared.joint functionality. Not useful.

plotSpec  = false;
rGroups   = [1 2];
extraOpts = assignopts(who, varargin);

yDims = params.yDims;
numGroups = length(yDims);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
  
obs_block_idxs = get_block_idxs(yDims); % Index rows of C
lat_block_idxs = get_block_idxs(xDim_across + xDim_within); % Index cols of C

% Reformat the loading matrix C to consolidate across-group dimensions
% into the same (appropiate) columns, then compute d_shared at the end.
C_across = []; 
d_shared.total = nan(1,numGroups);
shared_var.total = cell(1,numGroups);
cumulative_var.total = cell(1,numGroups);
for groupIdx = 1:numGroups
    
    % Indices of the current group observations
    currObsBlock = obs_block_idxs{groupIdx};
    obsIdxs = currObsBlock(1):currObsBlock(2);
    
    % Indices of current group latents
    currLatBlock = lat_block_idxs{groupIdx};
    latIdxs = currLatBlock(1):currLatBlock(2);
    
    % Indices of current across-group latents
    acrossLatIdxs = latIdxs(1:xDim_across);
    % Collect the specified across-group loading
    % dimensions for the current group
    C_across = [C_across; params.C(obsIdxs,acrossLatIdxs)];

    % Collect all latent dimensions for the current group
    C_total = params.C(obsIdxs,latIdxs);

    % Get the eigenvalue spectrum of the current group's loading matrix
    [~, D] = eig(C_total*C_total.');
    [spectrum,~] = sort(diag(D), 'descend');

    % Now find the smallest number of dimensions required to explain
    % shared variance within a certain cutoff
    shared_var.total{groupIdx} = spectrum ./ sum(spectrum);
    cumulative_var.total{groupIdx} = cumsum(spectrum) ./ sum(spectrum);
    [~, d_shared.total(groupIdx)] = max(cumulative_var.total{groupIdx} >= cutoffPC);

end

% Get the eigenvalue spectrum of the joint loading matrix
CC = C_across*C_across';

% Take the appropriate block of the shared covariance matrix.
obsBlock1 = obs_block_idxs{rGroups(1)};
obsIdxs1 = obsBlock1(1):obsBlock1(2);
obsBlock2 = obs_block_idxs{rGroups(2)};
obsIdxs2 = obsBlock2(1):obsBlock2(2);
crosscov = CC(obsIdxs1,obsIdxs2);

% Get the singular values of the across-group covariance matrix
s = svd(crosscov);
[s,~] = sort(s, 'descend'); % Just ensure that the singular values are sorted

% Now find the smallest number of dimensions required to explain
% shared across-group-covariance within a certain cutoff
shared_var.across = s ./ sum(s);
cumulative_var.across = cumsum(s) ./ sum(s);
[~, d_shared.across] = max(cumulative_var.across >= cutoffPC);

if plotSpec
    colors = generateColors(); % Generate custom plotting colors
    
    % Plot the eigenspectra and cumulative shared variance explained
    figure;
    numCol = 2;
    numRow = numGroups + 1;
    
    % Exclusivley across-group covariance
    subplot(numRow,numCol,1);
    hold on;
    plot(1:length(shared_var.across), shared_var.across, ...
         'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
    xlabel('Dimension, across-group');
    ylabel('Frac. shared cross-cov. exp.');
    xlim([0 length(shared_var.across)+1]);
    
    subplot(numRow,numCol,2);
    hold on;
    plot(1:length(cumulative_var.across), cumulative_var.across, ...
         'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
    line([0 length(cumulative_var.across)+1], [cutoffPC cutoffPC], ...
         'linestyle', '--', 'color', colors.reds{4}, 'linewidth', 1.5);
    plot(d_shared.across, cumulative_var.across(d_shared.across), ...
         'p', 'color', colors.reds{4}, 'markerfacecolor', ...
         colors.reds{4}, 'markersize', 10);
    xlim([0 length(cumulative_var.across)+1]);
    xlabel('No. dimensions, across-group');
    ylabel('Cumulative shared cross-cov. exp.');
    
    % For each group individually
    for groupIdx = 1:numGroups
        subplot(numRow,numCol,numCol*groupIdx+1);
        hold on;
        plot(1:length(shared_var.total{groupIdx}), shared_var.total{groupIdx}, ...
             'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
        xlabel(sprintf('Dimension, group %d', groupIdx));
        ylabel('Frac. shared var. exp.');
        xlim([0 length(shared_var.total{groupIdx})+1]);

        subplot(numRow,numCol,numCol*groupIdx+2);
        hold on;
        plot(1:length(cumulative_var.total{groupIdx}), cumulative_var.total{groupIdx}, ...
             'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
        line([0 length(cumulative_var.total{groupIdx})+1], [cutoffPC cutoffPC], ...
         'linestyle', '--', 'color', colors.reds{4}, 'linewidth', 1.5);
        plot(d_shared.total(groupIdx), cumulative_var.total{groupIdx}(d_shared.total(groupIdx)), ...
         'p', 'color', colors.reds{4}, 'markerfacecolor', colors.reds{4}, ...
         'markersize', 10);
        xlim([0 length(cumulative_var.total{groupIdx})+1]);
        xlabel(sprintf('No. dimensions, group %d', groupIdx));
        ylabel('Cumulative shared var. exp.');
        
    end
   
end
