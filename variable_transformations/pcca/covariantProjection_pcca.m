function Xcov = covariantProjection_pcca(Ys, params, varargin)
%
% Xcov = covariantProjection_pcca(Ys, params, ...)
%
% Description: Project the latents inferred by a pCCA model onto
%              covariant modes, which capture the maximal covariance
%              between a pair of groups.
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
%                  analyze. Order does not change the computation, but it
%                  does change the order of groups in the outputs. 
%                  (default: [1 2])
%
% Outputs:
%
%     Xcov  -- (2*xDim x N) array; projections onto covariant modes. The
%              first xDim latents correspond to group groupIdxs(1). The
%              next xDim latents correspond to group groupIdxs(2).
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     26 Feb 2022 -- Initial full revision.

groupIdxs = [1 2];
assignopts(who,varargin);
numGroups = length(groupIdxs);

% Infer latents
[X, ~] = pcca_estep(Ys, params);

% Compute predictive modes
[~, V] = covariantModes_pcca(params, 'groupIdxs', groupIdxs);

% Project latents onto covariant modes for each group
C = params.Cs(groupIdxs);
xDim = size(C{1},2);
Xcov = nan(numGroups*xDim,size(X.mean,2));
for groupIdx = 1:numGroups
    % Project onto orthogonal basis
    Xcov((groupIdx-1)*xDim + (1:xDim),:) = V{groupIdx}' * C{groupIdx} * X.mean;
end
