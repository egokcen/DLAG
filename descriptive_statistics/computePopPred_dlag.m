function pred = computePopPred_dlag(params, varargin)
%
% pred = computePopPred_dlag(params, ...)
%
% Description: Compute:
%                - The 'population R^2' (a DLAG fit-derived coefficient of 
%                  determination) between a pair of groups.
%                - Predictive power within each predictive mode.
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
%                  analyze. Order matters: groupIdxs(1) gives the source
%                  group; groupIdxs(2) gives the target group.
%                  (default: [1 2])
%
% Outputs:
%
%     pred -- structure with the following fields:
%                 tpred      --  float; total predictive power
%                                tpred = sum(mpred).
%                 mpred      -- (1 x xDim_across) array; predictive power
%                               within each predictive mode.
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
pred.tpred = 0;
pred.mpred = 0;

% Only proceed if there are across-group dimensions
if xDim_across > 0
    
    % Initialize the rest of the output structure
    pred.mpred = nan(1,xDim_across);

    % Compute predictive modes
    [P, targetVar, ~, ~, ~] = predictiveModes_dlag(params,'groupIdxs',groupIdxs);

    % Predictive power (R^2) within each predictive mode
    pred.mpred = (diag(P)').^2 ./ targetVar;

    % Total 'population R^2'
    pred.tpred = sum(pred.mpred);

end