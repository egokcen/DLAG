function [P, targetVar, U, V, H] = predictiveModes_dlag(params, varargin)
%
% [P, targetVar, U, V, H] = predictiveModes_dlag(params, ...)
%
% Description: Compute the predictive modes between a source and target
%              group, which capture maximal 0-lag predictive power from
%              source to target.
%                  P = V{1}'*Sig_11^(-0.5)*Ca{1}*Ka*Ca{2}'*V{2}
%                    = U{1}'*Ca{1}*Ka*Ca{2}'*V{2}
%                    = H{1}*H{2}'
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
%     P  -- (xDim_across x xDim_across) array; diagonal matrix with the
%           cross-area predictive power along each predictive mode.
%     targetVar -- float; total variance of the target group.
%     U  -- (sourceDim x xDim_across) array;
%           predictive modes for the source group, which are uncorrelated 
%           but not necessarily orthogonal. U = Sig_11^(-.5)*V{1}
%     V  -- (1 x 2) cell array; V{i} -- (yDims(groupIdxs(i)) x xDim_across)
%           array; predictive modes for group groupIdxs(i), which are 
%           orthogonal but not necessarily uncorrelated.
%     H  -- (1 x 2) cell array; H{i} -- (yDims(groupIdxs(i)) x xDim_across)
%           array; Projection of latents onto predictive modes: 
%           H{1} = U'Ca{1}Ka; H{2} = V{2}'Ca{2}Ka
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     17 Mar 2021 -- Initial full revision.
%     20 Oct 2021 -- Updated documentation to clarify group order of
%                    outputs.

groupIdxs = [1 2];
assignopts(who,varargin);
numGroups = length(groupIdxs);

% Initialize outputs
P = [];
targetVar = nan;
U = [];
V = cell(1,numGroups);
H = cell(1,numGroups);

% Only proceed if there are across-group dimensions
xDim_across = params.xDim_across;
if xDim_across > 0

    % Get the 0-lag portion of the across-group GP kernel matrix, Ka(t,t)
    Ka = getZeroLagK_across(params,'groupIdxs',groupIdxs);

    % Get observation parameters for each group
    groupParams = partitionParams_dlag(params);
    % Take only the group pair of interest
    groupParams = groupParams(groupIdxs);
    xDim_within = params.xDim_within(groupIdxs);
    Ca = cell(1,numGroups);
    Cw = cell(1,numGroups);
    R = cell(1,numGroups);
    for groupIdx = 1:numGroups
        R{groupIdx} = groupParams{groupIdx}.R;
        acrossParams = getSubsetParams_dlag(groupParams{groupIdx}, 1:xDim_across, cell(1,numGroups));
        Ca{groupIdx} = acrossParams.C; 
        withinParams = getSubsetParams_dlag(groupParams{groupIdx},[], {1:xDim_within(groupIdx)});
        Cw{groupIdx} = withinParams.C;
    end

    % Compute least-squares matrix
    Sig_11 = Ca{1}*Ca{1}' + Cw{1}*Cw{1}' + R{1};
    Sig_22 = Ca{2}*Ca{2}' + Cw{2}*Cw{2}' + R{2};
    Sig_12 = Ca{1}*Ka*Ca{2}';
    Corr_12 = Sig_11^(-0.5) * Sig_12;

    % Compute predictive modes
    [V1, P, V2] = svd(Corr_12, 'econ');
    P = P(1:xDim_across,1:xDim_across);
    targetVar = trace(Sig_22);
    V{1} = V1(:,1:xDim_across);
    V{2} = V2(:,1:xDim_across);
    
    U = Sig_11^(-0.5)*V{1};
    
    H{1} = U'*Ca{1}*Ka^(.5);
    H{2} = V{2}'*Ca{2}*Ka^(.5);
    
end