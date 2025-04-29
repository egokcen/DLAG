function outparams = getSubsetParams_dlag(inparams, xDims_across, xDims_within)
%
% outparams = getSubsetParams_dlag(inparams, xDims_across, xDims_within)
%
% Description: Return a DLAG model with only parameters corresponding to
%              the specified within- and across-group parameters.
%
% Arguments:
%
%     inparams  -- Structure containing DLAG model parameters.
%                  Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
%                    eps_across   -- (1 x xDim_across) GP noise variances
%                    gamma_within -- (1 x numGroups) cell array; 
%                                    GP timescales for each group
%                    eps_within   -- (1 x numGroups) cell array;
%                                    GP noise variances for each group
%                    if covType == 'sg'
%                        nu_across -- (1 x xDim_across) array; center
%                                     frequencies for spectral Gaussians;
%                                     convert to 1/time via 
%                                     nu_across./binWidth 
%                        nu_within -- (1 x numGroups) cell array; 
%                                     center frequencies for each group
%                    d            -- (yDim x 1) array; observation mean
%                    C            -- (yDim x (numGroups*xDim)) array;
%                                    mapping between low- and high-d spaces
%                    R            -- (yDim x yDim) array; observation noise
%                                    covariance 
%                    DelayMatrix  -- (numGroups x xDim_across) array;
%                                    delays from across-group latents to 
%                                    observed variables. NOTE: Delays are
%                                    reported as (real-valued) number of
%                                    time-steps.
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%     xDims_across -- (1 x numDims) array; across-group state 
%                     dimensions to be retained in outparams.
%     xDims_within -- (1 x numGroups) cell array; each element is a vector 
%                     of within-group state dimensions to be retained in 
%                     outparams.
%
% Outputs:
%
%     outparams -- Structure containing subset of DLAG model parameters.
%                  Same format as inparams, above.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     14 May 2020 -- Initial full revision.
%     22 Feb 2023 -- Added spectral Gaussian compatibility.

yDims = inparams.yDims;
numGroups = length(yDims);

% Initialize output structure
outparams = inparams;
outparams.xDim_across = length(xDims_across);
for groupIdx = 1:numGroups
    outparams.xDim_within(groupIdx) = length(xDims_within{groupIdx}); 
end

% State parameters
outparams.gamma_across = outparams.gamma_across(xDims_across);
outparams.eps_across = outparams.eps_across(xDims_across);
outparams.DelayMatrix = outparams.DelayMatrix(:,xDims_across);
if isequal(outparams.covType,'sg')
    outparams.nu_across = outparams.nu_across(xDims_across); 
end
for groupIdx = 1:numGroups
    outparams.gamma_within{groupIdx} = outparams.gamma_within{groupIdx}(xDims_within{groupIdx});
    outparams.eps_within{groupIdx} = outparams.eps_within{groupIdx}(xDims_within{groupIdx});
    if isequal(outparams.covType,'sg')
        outparams.nu_within{groupIdx} = outparams.nu_within{groupIdx}(xDims_within{groupIdx}); 
    end
end

% Loading matrix
lat_block_idxs = get_block_idxs(inparams.xDim_across + inparams.xDim_within); % Index cols of C
keptIdxs = []; % Track which columns of C to keep
for groupIdx = 1:numGroups 

    % Indices of the current group's latents
    latBlock = lat_block_idxs{groupIdx};
    latIdxs = latBlock(1):latBlock(2);
    
    % Determine which columns to keep for this group
    keptIdxs = [keptIdxs latIdxs([xDims_across inparams.xDim_across + xDims_within{groupIdx}])];

end

% Retain only desired dimensions.
outparams.C = outparams.C(:,keptIdxs);