function Xcorr = correlativeProjection_pcca(Ys, params, varargin)
%
% Xcorr = correlativeProjection_pcca(Ys, params, ...)
%
% Description: Project the latents inferred by a pCCA model onto
%              correlative modes, which capture the maximal correlation 
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
%     orth      -- logical; If true, project onto an orthogonal (as opposed
%                  to uncorrelated) basis. (default: false)
%
% Outputs:
%
%     Xcorr -- (2*xDim x N) array; projections onto correlative modes. The
%              first xDim latents correspond to group groupIdxs(1). The
%              next xDim latents correspond to group groupIdxs(2).
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     20 Mar 2021 -- Initial full revision.
%     20 Oct 2021 -- Updated documentation to clarify group order in Xcorr.

groupIdxs = [1 2];
orth = false;
assignopts(who,varargin);
numGroups = length(groupIdxs);

% Infer latents
[X, ~] = pcca_estep(Ys, params);

% Compute predictive modes
[~, U, V] = correlativeModes_pcca(params, 'groupIdxs', groupIdxs);

% Project latents onto correlative modes for each group
C = params.Cs(groupIdxs);
xDim = size(C{1},2);
Xcorr = nan(numGroups*xDim,size(X.mean,2));
for groupIdx = 1:numGroups
    if orth
        % Project onto orthogonal basis
        Xcorr((groupIdx-1)*xDim + (1:xDim),:) = V{groupIdx}' * C{groupIdx} * X.mean;
    else
        % Project onto uncorrelated basis
        Xcorr((groupIdx-1)*xDim + (1:xDim),:) = U{groupIdx}' * C{groupIdx} * X.mean;
    end
end
