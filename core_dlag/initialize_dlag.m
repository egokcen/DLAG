function startParams = initialize_dlag(seqTrain,init_method,varargin)
%
% initialize_dlag
%
% Description: Initialize DLAG parameters for the EM algorithm.
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
%     init_method  -- string; Specify how DLAG parameters should be
%                     initialized. Options currently supported:
%                         'pCCA'   -- parameters are initialized using 
%                                     probabilistic canonical correlation 
%                                     analysis (pCCA), described in 
%                                     (Bach & Jordan, 2005).
%                         'params' -- parameters are initialized to a 
%                                     particular realization of DLAG
%                                     model parameters. For example, one 
%                                     may wish to continue training
%                                     a DLAG model for which training was 
%                                     cut short, or performed in a 
%                                     different context.
%
%     Optional:
%
%     xDims_across -- (1 x numDims) array; across-group state 
%                     dimensionalities to be modeled. For example, to fit
%                     3 separate models, each with 1, 2, and 3 across-group
%                     dimensions, input 1:3 or [1 2 3]. (default: [3])
%     xDims_within -- (1 x numGroups) cell array; each element is a vector of
%                     within-group state dimensionalities to be modeled, of
%                     the same format as xDims_across. For example, {1:3,
%                     1:3} specifies desired within-group state
%                     dimensionalities for data with two groups.
%                     If xDims_within contains only one element, but there
%                     are multiple groups (as indicated by yDims,
%                     xDims_across), then all groups will be assigned the
%                     same within-group dimensionalities specified in that
%                     vector. (default: {[1]})
%     yDims        -- (1 x numGroups) array; Specify the number features 
%                     (neurons) in each group (area). Elements in yDims
%                     should match the format of data in seqTrain and seqTest. 
%     rGroups      -- (1 x 2) array; Used to assess cross-validated
%                     performance via pairwise regression. Each element
%                     specifies a group to be included in the regression.
%                     (default: [1 2])
%     initParams   -- Structure containing DLAG model parameters. Contains
%                     same fields as 'startParams', documented below.
%                     Only relevant if 'init_params' is 'params'.
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
%                    xDim_within  -- (1 x numGroups) array; number 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.   

xDim_across   = 3;
xDim_within   = [];
yDims         = [];
rGroups       = [1 2];
initParams    = {};
parallelize   = false;
extraOpts     = assignopts(who, varargin);

if ~parallelize
    fprintf('Initializing DLAG using %s\n',init_method);
end

switch(init_method)
    case 'pCCA'
        startParams = init_pCCA_dlag(seqTrain, extraOpts{:}, ...
            'xDim_across', xDim_across, 'xDim_within', xDim_within, ...
            'yDims', yDims); 
    case 'params'
        startParams = initParams;
    otherwise
        fprintf('\nError: Invalid Option\n');
end
