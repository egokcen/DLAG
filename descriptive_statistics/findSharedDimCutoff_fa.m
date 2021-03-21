function d_shared = findSharedDimCutoff_fa(params, cutoffPC, varargin)
%
% d_shared = findSharedDimCutoff_fa(params, cutoffPC, ...)
%
% Description: Find the number of latent dimensions that explain a certain
%              percentage of shared variance.
%
% Arguments:
% 
%     Required:
%
%     params   -- Structure containing FA model parameters.
%                 Contains the (relevant) fields
% 
%                 C -- (yDim x xDim) array; factor loadings
%
%     cutoffPC -- float; cutoff proportion of shared variance explained 
%                 (between 0 and 1)
% 
%     Optional:
%
%     plotSpec -- logical; if true, plot the spectra and cumulative shared
%                 variance explained. (default: false)
%
% Outputs:
%
%     d_shared -- int; number of dimensions required to explain cutoffPC of
%                 the shared variance.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

plotSpec  = false;
extraOpts = assignopts(who, varargin);

% Get the eigenvalue spectrum of the shared covariance matrix
[~, D] = eig(params.C * params.C'); % (yDim x yDim) array
[spectrum,~] = sort(diag(D), 'descend');

% Now find the smallest number of dimensions required to explain
% shared variance within a certain cutoff
shared_var = spectrum ./ sum(spectrum);
cumulative_var = cumsum(spectrum) ./ sum(spectrum);
[~, d_shared] = max(cumulative_var >= cutoffPC);

if plotSpec
    colors = generateColors(); % Generate custom plotting colors
    
    % Plot the spectrum and cumulative shared across-group-covariance explained
    figure;
    subplot(1,2,1);
    hold on;
    plot(1:length(shared_var), shared_var, ...
         'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
    xlabel('Dimension');
    ylabel('Frac. shared var. exp.');
    xlim([0 length(shared_var)+1]);
    
    subplot(1,2,2);
    hold on;
    plot(1:length(cumulative_var), cumulative_var, ...
         'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
    line([0 length(cumulative_var)+1], [cutoffPC cutoffPC], ...
         'linestyle', '--', 'color', colors.reds{4}, 'linewidth', 1.5);
    plot(d_shared, cumulative_var(d_shared), ...
         'p', 'color', colors.reds{4}, 'markerfacecolor', ...
         colors.reds{4}, 'markersize', 10);
    xlim([0 length(cumulative_var)+1]);
    xlabel('No. dimensions');
    ylabel('Cumulative shared var. exp.');
   
end
