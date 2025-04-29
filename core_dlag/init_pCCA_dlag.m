function startParams = init_pCCA_dlag(seqTrain,varargin)
%
% startParams = init_pCCA_dlag(seqTrain,...)
%
% Description: Initialize DLAG parameters using parameters estimated 
%              by pCCA, setting delays to startDelay and timescales to 
%              startTau. Within-group dimensions are initialized as the 
%              dimensions that explain the most variance in each group, 
%              while also being uncorrelated to the dimensions returned by
%              pCCA.
%
% Arguments:
%
%     Required: 
%
%     seqTrain -- training data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
%
%     Optional:
%
%     xDim_across  -- int; number of across-group latent variables
%     xDim_within  -- (1 x numGroups) array; number of within-group latents
%                     in each group
%     yDims        -- (1 x numGroups) array; Specify the number features 
%                     (neurons) in each group (area). Elements in yDims
%                     should match the format of data in seqTrain and seqTest. 
%     binWidth     -- float; bin width or sample period, in units of time.
%                     (default: 20)
%                     Note: For all other temporal variables, keep units
%                     consistent with binWidth. 
%     startTau     -- float; Initial GP timescale, in units oftime 
%                     (default: 2*binWidth)
%     startEps     -- float; Initial GP noise variance (default: 1e-3)
%     startDelay   -- (1 x numGroups) array; Initial delays between 
%                     across-area latents and groups. All across-area 
%                     latents reach a particular group after the same
%                     specified delay. Different groups can have different
%                     delays. Entries in units of time. (default: [])
%     startNu      -- float; Initial GP center frequency (for spectral
%                     Gaussian kernels only), in units of 1/time (same
%                     time units as binWidth) (default: 1/(10*binWidth))
%     covType      -- string; Specify GP covariance kernel type. Options
%                     currently supported:
%                         'rbf' -- Radial basis function, or squared
%                                  exponential kernel
%                         'sg'  -- spectral Gaussian, or Gauss-cosine
%                                  kernel
%     parallelize  -- logical; Set to true to use Matlab's parfor construct
%                     to parallelize each fold and latent dimensionality 
%                     using multiple cores. (default: false)
%
% Outputs:
%
%     startParams -- Structure containing DLAG model parameters at which EM
%                    algorithm is initialized. Contains the fields
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
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.   
%     17 Apr 2020 -- Added 0-within-group dimension functionality
%     15 May 2020 -- Convert data to pCCA-compatible format with seq2cell2D
%     19 Feb 2023 -- Added spectral Gaussian compatibility

xDim_across   = 3;
xDim_within   = [];
yDims         = [];
binWidth      = 20;         % in units of time
startTau      = 2*binWidth; % in units of time
startDelay    = [];         % in units of time
startEps      = 1e-3;
startNu       = 1/(10*binWidth);
covType       = 'rbf';
parallelize   = false;
extraOpts     = assignopts(who, varargin);
numGroups     = length(yDims);

% ========================================
% Initialize observation model parameters
% ========================================

% Convert data into format compatible with pCCA
Ys = seq2cell2D(seqTrain, yDims, 'datafield', 'y');
% Fit pCCA model
[pccaParams, pccaLL] = em_pcca(Ys, xDim_across, ...
                               'parallelize', parallelize, extraOpts{:});

% Compute the within-group dimensions: Explain the most covariance in each
% area while remaining uncorrelatyDed to the pCCA dimensions
C = cell(1,numGroups); % Includes within- and across-group latents for each area
for groupIdx = 1:numGroups
    currY = Ys{groupIdx};
    C_across = pccaParams.Cs{groupIdx};
    C_uncorr = [];
    % C should only have within-group columns if xDim_within is non-zero
    if xDim_within(groupIdx) > 0
        % Compute the uncorrelated subspace
        covY = cov(currY');
        [~, ~, C_uncorr] = svd(C_across' * covY);
        C_uncorr = C_uncorr(:,xDim_across+1:xDim_across+xDim_within(groupIdx));
    end
    % Construct big C for this group
    C{groupIdx} = [C_across C_uncorr];
end
startParams.d = cat(1,pccaParams.ds{:});
% Take only the diagonal elements of pCCA's R matrix (which is
% block-diagonal)
startParams.R = diag(diag(blkdiag(pccaParams.Rs{:})));

% Restructure C matrix
startParams.C = blkdiag(C{:});

% Delay matrix
if isempty(startDelay)
    % No initial delays were specified. Initialize to zeros.
    startParams.DelayMatrix = zeros(numGroups,xDim_across);
else
    % All across-area latents reach a particular group after the same
    % specified delay. Different groups can have different delays.
    assert(length(startDelay) == numGroups);
    for groupIdx = 1:numGroups
        startParams.DelayMatrix(groupIdx,:) = ones(1,xDim_across) .* (startDelay(groupIdx) / binWidth);
    end
end
startParams.xDim_across = xDim_across;
startParams.xDim_within = xDim_within;
startParams.yDims = yDims;

% ==================================
% Initialize state model parameters
% ==================================
startParams.covType = covType;
% GP timescale
% Assume binWidth is the time step size.
startParams.gamma_across = (binWidth / startTau)^2 * ones(1, xDim_across);
startParams.gamma_within = cell(1,numGroups);
for groupIdx = 1:numGroups
    % Leave gamma_within empty if xDim_within is 0
    if xDim_within(groupIdx) > 0
        startParams.gamma_within{groupIdx} = (binWidth / startTau)^2 * ones(1, xDim_within(groupIdx));
    end
end

% GP noise variance
startParams.eps_across = startEps * ones(1, xDim_across);
startParams.eps_within = cell(1,numGroups);
for groupIdx = 1:numGroups
    % Leave eps_within empty if xDim_within is 0
    if xDim_within(groupIdx) > 0
        startParams.eps_within{groupIdx} = startEps * ones(1, xDim_within(groupIdx));
    end
end

% GP center frequency
if isequal(covType, 'sg')
    startParams.nu_across = startNu * binWidth * ones(1, xDim_across);
    startParams.nu_within = cell(1,numGroups);
    for groupIdx = 1:numGroups
        % Leave nu_within empty if xDim_within is 0
        if xDim_within(groupIdx) > 0
            startParams.nu_within{groupIdx} = startNu * binWidth * ones(1, xDim_within(groupIdx));
        end
    end
end

% Define parameter constraints
startParams.notes.learnKernelParams      = true;
startParams.notes.learnGPNoise           = false;
startParams.notes.RforceDiagonal         = true;
