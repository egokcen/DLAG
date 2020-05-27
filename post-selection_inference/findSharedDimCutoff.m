function d_shared = findSharedDimCutoff(params, cutoffPC)
%
% d_shared = findSharedDimCutoff(params, cutoffPC)
%
% Description: Find the number of across-group dimensions that explain
%              a certain percentage of across-group shared variance. Also
%              find the number of within-group dimensions that explain a
%              certain percentage of within-group shared variance, for each
%              group individually.
%
% Arguments:
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
% Outputs:
%
%     d_shared -- structure with the following fields:
%                 across -- int; number of dimensions required
%                           to explain cutoffPC of the across-group shared
%                           variance
%                 within -- (1 x numGroups) array; number of dimensions
%                           required to explain cutoffPC of the
%                           within-group shared variance for each group.
%                           within(i) will be NaN wherever xDim_within(i)
%                           is 0.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     11 Apr 2020 -- Initial full revision.
%     17 Apr 2020 -- Added 0-within-group dimension functionality
%     24 May 2020 -- Updated argument documentation.

yDims = params.yDims;
numGroups = length(yDims);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
  
obs_block_idxs = get_block_idxs(yDims); % Index rows of C
lat_block_idxs = get_block_idxs(xDim_across + xDim_within); % Index cols of C

% Stack across-group dimensions into one matrix, 
% and find d_shared at the end
C_across = []; 
d_shared.within = nan(1,numGroups);
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
    
    % Only deal with within-group latents if there are any.
    if xDim_within(groupIdx) > 0
        % Indices of current within-group latents
        withinLatIdxs = latIdxs(xDim_across+1:end);
        % Collect the specified within-group loading
        % dimensions for the current group
        C_within = params.C(obsIdxs,withinLatIdxs);

        % Get the eigenvalue spectrum of the current within-group loading matrix
        [~, D] = eig(C_within*C_within.');
        [spectrum,~] = sort(diag(D), 'descend');

        % Now find the smallest number of dimensions required to explain
        % within-group shared variance within a certain cutoff
        for d = 1:xDim_within(groupIdx)
            shared_var = sum(spectrum(1:d))/sum(spectrum);
            if shared_var >= cutoffPC
                break; 
            end
        end
        d_shared.within(groupIdx) = d;
    end
end

% Get the eigenvalue spectrum of the across-group loading matrix
[~, D] = eig(C_across*C_across.');
[spectrum,~] = sort(diag(D), 'descend');

% Now find the smallest number of dimensions required to explain
% across-group shared variance within a certain cutoff
for d = 1:xDim_across
    shared_var = sum(spectrum(1:d))/sum(spectrum);
    if shared_var >= cutoffPC
        break; 
    end
end
d_shared.across = d;
