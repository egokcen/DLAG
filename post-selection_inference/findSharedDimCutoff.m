function d_shared = findSharedDimCutoff(params, cutoffPC, varargin)
%
% d_shared = findSharedDimCutoff(params, cutoffPC, ...)
%
% Description: Find the number of overall dimensions (within or across) 
%              that explain a certain percentage of shared variance, 
%              jointly across all groups and for each group individually.
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
%
% Outputs:
%
%     d_shared -- structure with the following fields:
%                 joint -- int; number of dimensions required
%                           to explain cutoffPC of the joint shared
%                           variance across all groups
%                 indiv -- (1 x numGroups) array; number of dimensions
%                           required to explain cutoffPC of the
%                           shared variance for each group individually.
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

plotSpec = false;
extraOpts     = assignopts(who, varargin);

yDims = params.yDims;
numGroups = length(yDims);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
  
obs_block_idxs = get_block_idxs(yDims); % Index rows of C
lat_block_idxs = get_block_idxs(xDim_across + xDim_within); % Index cols of C

% Reformat the loading matrix C to consolidate across-group dimensions
% into the same (appropiate) columns, then compute d_shared at the end.
C_across = []; 
d_shared.indiv = nan(1,numGroups);
shared_var.indiv = cell(1,numGroups);
cumulative_var.indiv = cell(1,numGroups);
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
    C_indiv = params.C(obsIdxs,latIdxs);

    % Get the eigenvalue spectrum of the current group's loading matrix
    [~, D] = eig(C_indiv*C_indiv.');
    [spectrum,~] = sort(diag(D), 'descend');

    % Now find the smallest number of dimensions required to explain
    % shared variance within a certain cutoff
    shared_var.indiv{groupIdx} = spectrum ./ sum(spectrum);
    cumulative_var.indiv{groupIdx} = cumsum(spectrum) ./ sum(spectrum);
    [~, d_shared.indiv(groupIdx)] = max(cumulative_var.indiv{groupIdx} >= cutoffPC);

end

% Construct the restructured joint loading matrix
C_joint = C_across;
for groupIdx = 1:numGroups
    if xDim_within(groupIdx) > 0
        % Indices of current group latents
        currLatBlock = lat_block_idxs{groupIdx};
        latIdxs = currLatBlock(1):currLatBlock(2);

        % Indices of current within-group latents
        withinLatIdxs = latIdxs(xDim_across+1:end);
        
        C_joint = [C_joint params.C(:,withinLatIdxs)];
    end
end

% Get the eigenvalue spectrum of the across-group loading matrix
[~, D] = eig(C_joint*C_joint.');
[spectrum,~] = sort(diag(D), 'descend');

% Now find the smallest number of dimensions required to explain
% shared variance within a certain cutoff
shared_var.joint = spectrum ./ sum(spectrum);
cumulative_var.joint = cumsum(spectrum) ./ sum(spectrum);
[~, d_shared.joint] = max(cumulative_var.joint >= cutoffPC);

if plotSpec
    colors = generateColors(); % Generate custom plotting colors
    
    % Plot the eigenspectra and cumulative shared variance explained
    figure;
    numCol = 2;
    numRow = numGroups + 1;
    
    % Jointly across groups
    subplot(numRow,numCol,1);
    hold on;
    plot(1:length(shared_var.joint), shared_var.joint, ...
         'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
    xlabel('Dimension, all groups');
    ylabel('Frac. shared var. exp.');
    xlim([0 length(shared_var.joint)+1]);
    
    subplot(numRow,numCol,2);
    hold on;
    plot(1:length(cumulative_var.joint), cumulative_var.joint, ...
         'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
    line([0 length(cumulative_var.joint)+1], [cutoffPC cutoffPC], ...
         'linestyle', '--', 'color', colors.reds{4}, 'linewidth', 1.5);
    plot(d_shared.joint, cumulative_var.joint(d_shared.joint), ...
         'p', 'color', colors.reds{4}, 'markerfacecolor', ...
         colors.reds{4}, 'markersize', 10);
    xlim([0 length(cumulative_var.joint)+1]);
    xlabel('No. dimensions, all groups');
    ylabel('Cumulative shared var. exp.');
    
    % For each group individually
    for groupIdx = 1:numGroups
        subplot(numRow,numCol,numCol*groupIdx+1);
        hold on;
        plot(1:length(shared_var.indiv{groupIdx}), shared_var.indiv{groupIdx}, ...
             'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
        xlabel(sprintf('Dimension, group %d', groupIdx));
        ylabel('Frac. shared var. exp.');
        xlim([0 length(shared_var.indiv{groupIdx})+1]);

        subplot(numRow,numCol,numCol*groupIdx+2);
        hold on;
        plot(1:length(cumulative_var.indiv{groupIdx}), cumulative_var.indiv{groupIdx}, ...
             'o-', 'color', 'k', 'markerfacecolor', 'k', 'linewidth', 1.5);
        line([0 length(cumulative_var.indiv{groupIdx})+1], [cutoffPC cutoffPC], ...
         'linestyle', '--', 'color', colors.reds{4}, 'linewidth', 1.5);
        plot(d_shared.indiv(groupIdx), cumulative_var.indiv{groupIdx}(d_shared.indiv(groupIdx)), ...
         'p', 'color', colors.reds{4}, 'markerfacecolor', colors.reds{4}, ...
         'markersize', 10);
        xlim([0 length(cumulative_var.indiv{groupIdx})+1]);
        xlabel(sprintf('No. dimensions, group %d', groupIdx));
        ylabel('Cumulative shared var. exp.');
        
    end
   
end
