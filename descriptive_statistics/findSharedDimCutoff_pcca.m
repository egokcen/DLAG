function d_shared = findSharedDimCutoff_pcca(params, cutoffPC, varargin)
%
% d_shared = findSharedDimCutoff_pcca(params, cutoffPC, ...)
%
% Description: Find the number of across-group dimensions 
%              that explain a certain percentage of shared 
%              across-group-covariance.
%
% Arguments:
% 
%     Required:
%
%     params   -- Structure containing pCCA model parameters.
%                 Contains the (relevant) fields
% 
%                 Cs -- (1 x numGroups) cell array; List of factor loadings 
%                       {(y1Dim x xDim), (y2Dim x xDim), ...}
%
%     cutoffPC -- float; cutoff proportion of shared across-group-
%                 covariance explained (between 0 and 1)
% 
%     Optional:
%
%     plotSpec -- logical; if true, plot the spectra and cumulative shared
%                 across-group-covariance explained. (default: false)
%     rGroups  -- (1 x 2) int array; Indexes of the two groups whose cross-
%                 covariance we wish to consider. (default: [1 2])
%
% Outputs:
%
%     d_shared -- int; number of dimensions required to explain cutoffPC of
%                 the shared across-group-covariance.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

plotSpec  = false;
rGroups   = [1 2];
extraOpts = assignopts(who, varargin);

numGroups = length(params.Cs);
xDim = size(params.Cs{1},2);
yDims = nan(1,numGroups);
for groupIdx = 1:numGroups
    yDims(groupIdx) = size(params.Cs{groupIdx},1); 
end
  
obs_block_idxs = get_block_idxs(yDims); % Index rows of C

% Construct the full loading matrix
C = vertcat(params.Cs{:}); % (yDim x xDim) array

% Compute the shared covariance matrix
CC = C * C'; % (yDim x yDim) array

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
shared_cov = s ./ sum(s);
cumulative_cov = cumsum(s) ./ sum(s);
[~, d_shared] = max(cumulative_cov >= cutoffPC);

if plotSpec
    colors = generateColors(); % Generate custom plotting colors
    
    % Plot the spectrum and cumulative shared across-group-covariance explained
    figure;
    subplot(1,2,1);
    hold on;
    plot(1:length(shared_cov), shared_cov, ...
         'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
    xlabel('Dimension');
    ylabel('Frac. shared cross-cov. exp.');
    xlim([0 length(shared_cov)+1]);
    
    subplot(1,2,2);
    hold on;
    plot(1:length(cumulative_cov), cumulative_cov, ...
         'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
    line([0 length(cumulative_cov)+1], [cutoffPC cutoffPC], ...
         'linestyle', '--', 'color', colors.reds{4}, 'linewidth', 1.5);
    plot(d_shared, cumulative_cov(d_shared), ...
         'p', 'color', colors.reds{4}, 'markerfacecolor', ...
         colors.reds{4}, 'markersize', 10);
    xlim([0 length(cumulative_cov)+1]);
    xlabel('No. dimensions');
    ylabel('Cumulative shared cross-cov. exp.');
   
end
