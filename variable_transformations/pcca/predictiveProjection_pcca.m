function Xpred = predictiveProjection_pcca(Ys, params, varargin)
%
% Xpred = predictiveProjection_pcca(Ys, params, ...)
%
% Description: Project the latents inferred by a pCCA model onto predictive
%              modes, which capture the maximal predictive power between a 
%              source and target group.
%
% Arguments:
%
%     Required:
%
%     Ys      -- (1 x numGroups) cell array; list of data matrices 
%                {(y1Dim x N), (y2Dim x N), ...}
%     params  -- learned pCCA parameters (structure with fields Cs, Rs, ds)
%
%     Optional:
%
%     groupIdxs -- (1 x 2) int array; Specify which pair of groups to
%                  analyze. Order matters: groupIdxs(1) gives the source
%                  group; groupIdxs(2) gives the target group.
%                  (default: [1 2])
%     orth      -- logical; If true, project source latents onto an 
%                  orthogonal (as opposed) to uncorrelated) basis. 
%                  (default: false)
%
% Outputs:
%
%     Xpred -- (2*xDim x N) array; projections onto predictive modes
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Mar 2019 -- Initial full revision.
%     20 Mar 2021 -- Overhauled to use predictiveModes_pcca.m

groupIdxs = [1 2];
orth = false;
assignopts(who,varargin);
numGroups = length(groupIdxs);

% Infer latents
[X, ~] = pcca_estep(Ys, params);

% Compute predictive modes
[~, U, V] = predictiveModes_pcca(params, 'groupIdxs', groupIdxs);

% Project latents onto predictive modes for each group
C = params.Cs(groupIdxs);
xDim = size(C{1},2);
Xpred = nan(numGroups*xDim,size(X.mean,2));
for groupIdx = 1:numGroups
    if orth || groupIdx > 1
        % Project target and/or source onto orthogonal basis
        Xpred((groupIdx-1)*xDim + (1:xDim),:) = V{groupIdx}' * C{groupIdx} * X.mean;
    else
        % Project source onto uncorrelated basis
        Xpred((groupIdx-1)*xDim + (1:xDim),:) = U' * C{groupIdx} * X.mean;
    end
end
