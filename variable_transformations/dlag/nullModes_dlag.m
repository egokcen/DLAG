function [S, Uout, Vout] = nullModes_dlag(params, Uin)
%
% [S, Uout, Vout] = nullModes_dlag(params, Uin)
%
% Description: Compute the dominant modes of each group, after projecting
%              onto the null space of Uin.
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
%                    xDim_within  -- (1 x numGroups) array; number of 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
%     Uin -- (1 x 2) cell array; Uin{i} -- (yDims(i) x r) array; a basis 
%            for a subspace in group i.
%
% Outputs:
%
%     S    -- (1 x numGroups) cell array; S{i} -- (xDim-r x xDim-r) array; 
%             diagonal matrix with the singular values of each null mode.
%     Uout -- (1 x 2) cell array; U{i} -- (yDims(i) x xDim-r) array;
%             null modes (left singular vectors) for group i.
%     Vout -- (1 x 2) cell array; V{i} -- (yDims(i) x xDim-r) array;
%             right singular vectors for group i
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     22 Feb 2023 -- Initial full revision.

numGroups = length(params.yDims);

% Initialize outputs
S = cell(1,numGroups);
Uout = cell(1,numGroups);
Vout = cell(1,numGroups);

% Get observation parameters for each group
groupParams = partitionParams_dlag(params);

for groupIdx = 1:numGroups
    
    r = size(Uin{groupIdx},2);
    
    % Project onto input basis
    C = groupParams{groupIdx}.C;
    Qin = C' * Uin{groupIdx};      % xDim x r
    
    % Compute the null space in this xDim-dimensional space
    [Qnull,~] = qr(Qin);
    Qnull = Qnull(:,r+1:end);     % xDim x xDim-r
    
    % Remove Uin from the subspace spanned by C
    Cnull = C * (Qnull * Qnull');
    
    % Compute the dominant modes of the null space
    xDim = size(Cnull,2);
    [Ui, Si, Vi] = svd(Cnull, 'econ');
    S{groupIdx} = Si(1:(xDim-r),1:(xDim-r));
    Uout{groupIdx} = Ui(:,1:(xDim-r));
    Vout{groupIdx} = Vi(:,1:(xDim-r));
    
end
