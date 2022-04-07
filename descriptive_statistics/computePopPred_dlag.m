function pred = computePopPred_dlag(params, varargin)
%
% pred = computePopPred_dlag(params, ...)
%
% Description: Compute:
%                - The 'population R^2' (a DLAG fit-derived coefficient of 
%                  determination) between a pair of groups.
%                - Predictive power within each predictive mode.
%                - An aggregrate leave-group-out R^2 value.
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
%     zerolag   -- logical; set true to compute zero-lag modes, false
%                  to compute modes that factor in delays (default: true)
%
% Outputs:
%
%     pred -- structure with the following fields:
%                 tpred      -- (numGroups x 1) array; tpred(i) gives the
%                               total predictive power, sum(mpred(i,:)) 
%                               when group groupIdxs(i) is predicted by the
%                               other group.
%                 mpred      -- (numGroups x xDim_across) array; mpred(i,:)
%                               gives the predictive power within each 
%                               predictive mode, when group groupIdxs(i)
%                               is predicted by the other group (note that
%                               modes change for each group)
%                 R2lpo      -- float; R^2 value aggregated over all
%                               groups.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 Mar 2021 -- Initial full revision.
%     28 Feb 2022 -- Added aggregate leave-group-out metric.
%     06 Apr 2022 -- Added zerolag option.

groupIdxs = [1 2];
zerolag = true;
assignopts(who,varargin);

% Initialize part of the output structure
xDim_across = params.xDim_across;
numGroups = length(params.yDims);
pred.tpred = zeros(numGroups,1);
pred.mpred = zeros(numGroups,1);
pred.R2lpo = 0;

% Only proceed if there are across-group dimensions
if xDim_across > 0
    
    % Initialize the rest of the output structure
    pred.mpred = nan(numGroups,xDim_across);
    
    predVars = nan(numGroups,1);
    targetVars = nan(numGroups,1);

    for groupIdx = 1:length(groupIdxs)
        
        % Index of source group
        currSourceGroup = groupIdxs(groupIdx);
        
        % Index of the "held-out" target group
        currTargetGroup = setdiff(groupIdxs,currSourceGroup);
        
        % Compute predictive modes
        [P, targetVar, ~, ~, ~] = predictiveModes_dlag(params,'groupIdxs',[currSourceGroup currTargetGroup],'zerolag',zerolag);
        predVars(groupIdx) = trace(P.^2);
        targetVars(groupIdx) = targetVar;
        
        % Predictive power (R^2) within each predictive mode
        pred.mpred(groupIdx,:) = (diag(P)').^2 ./ targetVar;

        % Total 'population R^2'
        pred.tpred(groupIdx) = sum(pred.mpred(groupIdx,:));
    end
    
    % Aggregate R^2 across all groups
    pred.R2lpo = sum(predVars) / sum(targetVars);

end