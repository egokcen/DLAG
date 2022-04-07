function popcov = computePopCov_dlag(params, varargin)
%
% popcov = computePopCov_dlag(params, ...)
%
% Description: Compute:
%                - The total population covariance between a pair of groups.
%                - The covariance within each across-group latent
%                  dimension.
%                - The covariance within each covariant mode.
%                - The proportion of each covariant mode explained by
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
%     zerolag   -- logical; set true to compute zero-lag modes, false
%                  to compute modes that factor in delays (default: true)
%
% Outputs:
%
%     popcov -- structure with the following fields:
%                    tcov      --  float; total population covariance
%                                  tcorr = sum(mcov).
%                    lcov      -- (1 x xDim_across) array; covariance
%                                 within each across-group latent dimension
%                    mcov      -- (1 x xDim_across) array; covariance
%                                 within each covariant mode.
%                    mcovfrac  -- (xDim_across x xDim_across) array; 
%                                 mcovfrac(i,:) gives the proportion of 
%                                 mode i's covariance explained by each
%                                 across-group latent.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     26 Feb 2022 -- Initial full revision.
%     06 Apr 2022 -- Added zerolag option.

groupIdxs = [1 2];
zerolag = true;
assignopts(who,varargin);

% Initialize part of the output structure
xDim_across = params.xDim_across;
popcov.tcov = 0;

% Only proceed if there are across-group dimensions
if xDim_across > 0
    
    % Initialize the rest of the output structure
    popcov.lcov = nan(1,xDim_across);
    popcov.mcov = nan(1,xDim_across);
    popcov.mcovfrac = nan(xDim_across,xDim_across);

    % Compute covariant modes
    [P, ~, H] = covariantModes_dlag(params,'groupIdxs',groupIdxs,'zerolag',zerolag);

    % Covariance within each covariant mode
    popcov.mcov = diag(P)';

    % Total population covariance
    popcov.tcov = sum(popcov.mcov);

    for j = 1:xDim_across
        % Proportion of each mode explained by latent j
        popcov.mcovfrac(:,j) = H{1}(:,j).*H{2}(:,j) ./ popcov.mcov';
        % Covariance within latent j
        popcov.lcov(j) = sum(H{1}(:,j).*H{2}(:,j));
    end
    
else
    
    % Initialize the rest of the output structure
    popcov.lcov = 0;      % No covariance
    popcov.mcov = 0;      % No covariance
    popcov.mcovfrac = nan; % Fraction of covariance undefined

end