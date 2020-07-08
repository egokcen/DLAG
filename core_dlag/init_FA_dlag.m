function startParams = init_FA_dlag(seqTrain,varargin)
%
% startParams = init_FA_dlag(seqTrain,...)
%
% Description: Initialize DLAG parameters using parameters estimated 
%              by FA models fit independently to each group. Set timescales
%              to startTau. Across-group parameters will be empty.
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
%     covType      -- string; Specify GP covariance kernel type. Options
%                     currently supported:
%                         'rbf' -- Radial basis function, or squared
%                                  exponential kernel
%                     (default: 'rbf')
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
%     27 Jun 2020 -- Initial full revision.   

xDim_across   = 0;
xDim_within   = [];
yDims         = [];
binWidth      = 20;         % in units of time
startTau      = 2*binWidth; % in units of time
startDelay    = [];         % in units of time
startEps      = 1e-3;
covType       = 'rbf';
parallelize   = false;
extraOpts     = assignopts(who, varargin);
numGroups     = length(yDims);
xDim          = xDim_across + xDim_within;

% ========================================
% Initialize observation model parameters
% ========================================

% Convert data into format compatible with static methods
Ys = seq2pcca(seqTrain, yDims, 'datafield', 'y');
% Fit independent FA models to each group
Cs = cell(1,numGroups);
Rs = cell(1,numGroups);
ds = cell(1,numGroups);
for groupIdx = 1:numGroups
    fprintf('\nFitting FA model to group %d of %d\n', groupIdx, numGroups);
    [faParams, ~] = em_fa(Ys{groupIdx}, xDim(groupIdx), ...
                          'parallelize', parallelize, extraOpts{:});
    Cs{groupIdx} = faParams.C;
    Rs{groupIdx} = diag(faParams.R);
    ds{groupIdx} = faParams.d;
end
startParams.d = cat(1,ds{:});
startParams.R = blkdiag(Rs{:});
startParams.C = blkdiag(Cs{:});

% Delay matrix
startParams.DelayMatrix = [];
% Leave DelayMatrix empty if xDim_across is 0
if xDim_across > 0
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
% Leave gamma_across empty if xDim_across is 0
startParams.gamma_across = [];
if xDim_across > 0
    startParams.gamma_across = (binWidth / startTau)^2 * ones(1, xDim_across);
end
startParams.gamma_within = cell(1,numGroups);
for groupIdx = 1:numGroups
    % Leave gamma_within empty if xDim_within is 0
    if xDim_within(groupIdx) > 0
        startParams.gamma_within{groupIdx} = (binWidth / startTau)^2 * ones(1, xDim_within(groupIdx));
    end
end

% GP noise variance
% Leave eps_across empty if xDim_across is 0
startParams.eps_across = [];
if xDim_across > 0
    startParams.eps_across = startEps * ones(1, xDim_across);
end
startParams.eps_within = cell(1,numGroups);
for groupIdx = 1:numGroups
    % Leave eps_within empty if xDim_within is 0
    if xDim_within(groupIdx) > 0
        startParams.eps_within{groupIdx} = startEps * ones(1, xDim_within(groupIdx));
    end
end

% Define parameter constraints
startParams.notes.learnKernelParams      = true;
startParams.notes.learnGPNoise           = false;
startParams.notes.RforceDiagonal         = true;
