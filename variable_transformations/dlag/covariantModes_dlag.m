function [P, V, H] = covariantModes_dlag(params, varargin)
%
% [P, V, H] = covariantModes_dlag(params, ...)
%
% Description: Compute the covariant modes between a pair of groups, which
%              capture the maximal (0-lag) covariance across groups.
%                  P = V{1}'*Ca{1}*Ka*Ca{2}'*V{2}
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
%                  analyze. Order does not change the computation, but it
%                  does change the order of groups in the outputs. 
%                  (default: [1 2])
%     zerolag   -- logical; set true to compute zero-lag modes, false
%                  to compute modes that factor in delays (default: true)
%
% Outputs:
%
%     P -- (xDim_across x xDim_across) array; diagonal matrix with the
%          cross-area covariance along each covariant mode.
%     V -- (1 x 2) cell array; V{i} -- (yDims(groupIdxs(i)) x xDim_across)
%          array; covariant modes for group groupIdxs(i), which are 
%          orthogonal but not necessarily uncorrelated.
%     H -- (1 x 2) cell array; H{i} -- (yDims(groupIdxs(i)) x xDim_across)
%          array; Projection of latents onto covariant modes: 
%          H{i} = V{i}'Ca{i}Ka^(0.5)
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
numGroups = length(groupIdxs);

% Initialize outputs
P = [];
V = cell(1,numGroups);
H = cell(1,numGroups);

% Only proceed if there are across-group dimensions
xDim_across = params.xDim_across;
if xDim_across > 0

    if zerolag
        % Get the 0-lag portion of the across-group GP kernel matrix, Ka(t,t)
        Ka = getZeroLagK_across(params,'groupIdxs',groupIdxs);
    else
        % Otherwise, set Ka to the identity matrix, thereby maintaining
        % time delays in the computation of modes.
        Ka = eye(xDim_across);
    end

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

    % Compute cross-covariance matrix
    Sig_12 = Ca{1}*Ka*Ca{2}';

    % Compute covariant modes
    [V1, P, V2] = svd(Sig_12, 'econ');
    P = P(1:xDim_across,1:xDim_across);
    V{1} = V1(:,1:xDim_across);
    V{2} = V2(:,1:xDim_across);
    
    H{1} = V{1}'*Ca{1}*Ka^(.5);
    H{2} = V{2}'*Ca{2}*Ka^(.5);
    
end