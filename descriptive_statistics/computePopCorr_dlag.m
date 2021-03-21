function popcorr = computePopCorr_dlag(params, varargin)
%
% popcorr = computePopCorr_dlag(params, ...)
%
% Description: Compute:
%                - The total population correlation between a pair of groups.
%                - The correlation within each across-group latent
%                  dimension.
%                - The correlation within each correlative mode.
%                - The proportion of each correlative mode explained by
%                  each across-group latent.
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
%     groupIdxs -- (1 x 2) int array; Specify which pair of groups to
%                  analyze. Order doesn't matter. (default: [1 2])
%
% Outputs:
%
%     popcorr -- structure with the following fields:
%                    tcorr     --  float; total population correlation
%                                  tcorr = sum(mcorr).
%                    lcorr     -- (1 x xDim_across) array; correlation
%                                 within each across-group latent dimension
%                    mcorr     -- (1 x xDim_across) array; correlation
%                                 within each correlative mode.
%                    mcorrfrac -- (xDim_across x xDim_across) array; 
%                                 mcorrfrac(i,:) gives the proportion of 
%                                 mode i's correlation explained by each
%                                 across-group latent.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 Mar 2021 -- Initial full revision.

groupIdxs = [1 2];
assignopts(who,varargin);

% Initialize part of the output structure
xDim_across = params.xDim_across;
popcorr.tcorr = 0;

% Only proceed if there are across-group dimensions
if xDim_across > 0
    
    % Initialize the rest of the output structure
    popcorr.lcorr = nan(1,xDim_across);
    popcorr.mcorr = nan(1,xDim_across);
    popcorr.mcorrfrac = nan(xDim_across,xDim_across);

    % Compute correlative modes
    [P, ~, ~, H] = correlativeModes_dlag(params,'groupIdxs',groupIdxs);

    % Correlation within each correlative mode
    popcorr.mcorr = diag(P)';

    % Total population correlation
    popcorr.tcorr = sum(popcorr.mcorr);

    for j = 1:xDim_across
        % Proportion of each mode explained by latent j
        popcorr.mcorrfrac(:,j) = H{1}(:,j).*H{2}(:,j) ./ popcorr.mcorr';
        % Correlation within latent j
        popcorr.lcorr(j) = sum(H{1}(:,j).*H{2}(:,j));
    end
    
else
    
    % Initialize the rest of the output structure
    popcorr.lcorr = 0;      % No correlation
    popcorr.mcorr = 0;      % No correlation
    popcorr.mcorrfrac = nan; % Fraction of correlation undefined

end