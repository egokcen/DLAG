function [S, U, V, H] = dominantModes_dlag(params, varargin)
%
% [S, U, V, H] = dominantModes_dlag(params, ...)
%
% Description: Compute the dominant modes of each group, which
%              capture the greatest shared variance explained.
%                  C{i} = U{i}*S{i}*V{i}'
%
% Arguments:
%
%     Required:
%
%     params  -- Structure containing DLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
%                    eps_across   -- (1 x xDim_across) GP noise variances
%                    gamma_within -- (1 x numGroups) cell array; 
%                                    GP timescales for each group
%                    eps_within   -- (1 x numGroups) cell array;
%                                    GP noise variances for each group
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
%                    xDim_within  -- (1 x numGroups) array; number of 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
%     Optional:
%
%     includeAcross -- logical; set to true to include across-group latents
%                      in the computation (default: true)
%     includeWithin -- logical; set to true to include within-group latents
%                      in the computation (default: true)
%
% Outputs:
%
%     S -- (1 x numGroups) cell array; S{i} -- (xDim x xDim) array; 
%          diagonal matrix with the singular values of each dominant mode.
%     U -- (1 x 2) cell array; U{i} -- (yDims(i) x xDim) array;
%          dominant modes (left singular vectors) for group i.
%     V -- (1 x 2) cell array; V{i} -- (yDims(i) x xDim) array;
%          right singular vectors for group i
%     H -- (1 x 2) cell array; H{i} -- (yDims(i) x xDim_across) array;
%          Projection of latents onto dominant modes: 
%          H{i} = U{i}'*C{i}
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2021 -- Initial full revision.

includeAcross = true;
includeWithin = true;
assignopts(who,varargin);
assert(includeAcross || includeWithin);
numGroups = length(params.yDims);

% Initialize outputs
S = cell(1,numGroups);
U = cell(1,numGroups);
V = cell(1,numGroups);
H = cell(1,numGroups);

% Get observation parameters for each group
groupParams = partitionParams_dlag(params);

for groupIdx = 1:numGroups
    
    % Collect the appropriate within- and/or across-area loadings
    C = [];
    if includeAcross
        xDim_across = groupParams{groupIdx}.xDim_across;
        acrossParams = getSubsetParams_dlag(groupParams{groupIdx}, 1:xDim_across, cell(1,numGroups));
        C = [C acrossParams.C];
    end
    if includeWithin
        xDim_within = groupParams{groupIdx}.xDim_within;
        withinParams = getSubsetParams_dlag(groupParams{groupIdx},[], {1:xDim_within});
        C = [C withinParams.C];
    end
    
    % Compute dominant modes
    xDim = size(C,2);
    [Ui, Si, Vi] = svd(C, 'econ');
    S{groupIdx} = Si(1:xDim,1:xDim);
    U{groupIdx} = Ui(:,1:xDim);
    V{groupIdx} = Vi(:,1:xDim);
    H{groupIdx} = U{groupIdx}'*C;
    
end
