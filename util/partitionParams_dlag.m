function outparams = partitionParams_dlag(inparams)
%
% outparams = partitionParams_dlag(inparams)
%
% Description: Take a DLAG model with a given number of observation groups 
%              (numGroups) and return numGroups separate DLAG models, each
%              with one observation group.
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
%
% Outputs:
%
%     outparams -- (1 x numGroups) cell array; outparams(i) contains 
%                  the parameters of a DLAG model corresponding only to 
%                  observation group (i). Same format as inparams, above.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 May 2020 -- Initial full revision.
%     28 Jun 2020 -- Updated for 0-across-group dimension compatibility.
%     22 Feb 2023 -- Added spectral Gaussian compatibility.

numGroups = length(inparams.yDims);

% Indexes into observation and latent blocks, corresponding to groups
obs_block_idxs = get_block_idxs(inparams.yDims);
lat_block_idxs = get_block_idxs(inparams.xDim_across + inparams.xDim_within);

% Fill in output structure. Take DLAG parameters corresponding to single
% groups.
outparams = cell(1,numGroups);
for groupIdx = 1:numGroups
    outparams{groupIdx} = inparams;
    outparams{groupIdx}.xDim_across = inparams.xDim_across;
    outparams{groupIdx}.xDim_within = inparams.xDim_within(groupIdx); 
    outparams{groupIdx}.yDims = inparams.yDims(groupIdx);
    if inparams.xDim_across > 0
        outparams{groupIdx}.DelayMatrix = inparams.DelayMatrix(groupIdx,:);
    else
        outparams{groupIdx}.DelayMatrix = [];
    end
    outparams{groupIdx}.gamma_across = inparams.gamma_across;
    outparams{groupIdx}.eps_across = inparams.eps_across;
    outparams{groupIdx}.gamma_within = inparams.gamma_within(groupIdx);
    outparams{groupIdx}.eps_within = inparams.eps_within(groupIdx);
    if isequal(outparams{groupIdx}.covType,'sg')
        outparams{groupIdx}.nu_across = inparams.nu_across;
        outparams{groupIdx}.nu_within = inparams.nu_within(groupIdx);
    end
    
    % Observation model parameters require special indexing
    obsBlock = obs_block_idxs{groupIdx};
    obsIdxs = obsBlock(1):obsBlock(2);
    latBlock = lat_block_idxs{groupIdx};
    latIdxs = latBlock(1):latBlock(2);
    
    outparams{groupIdx}.d = inparams.d(obsIdxs);
    outparams{groupIdx}.R = inparams.R(obsIdxs,obsIdxs);
    outparams{groupIdx}.C = inparams.C(obsIdxs,latIdxs);
    
end