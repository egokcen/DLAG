function [P, V] = covariantModes_pcca(params, varargin)
%
% [P, V] = covariantModes_pcca(params, ...)
%
% Description: Compute the covariant modes between a pair of groups, 
%              which capture the maximal 0-lag covariance across groups.
%                  P = V{1}'*C{1}*C{2}'*V{2}
%
% Arguments:
%
%     Required:
%
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
%     P  -- (xDim x xDim) array; diagonal matrix with the
%           cross-area correlation along each correlative mode.
%     V  -- (1 x 2) cell array; V{i} -- (yDims(groupIdxs(i)) x xDim) array;
%           covariant modes for group groupIdxs(i), which are orthogonal
%           but not necessarily uncorrelated.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     26 Feb 2022 -- Initial full revision.

groupIdxs = [1 2];
assignopts(who,varargin);
numGroups = length(groupIdxs);

% Initialize outputs
V = cell(1,numGroups);

% Take only the parameters for the groups of interest
C = params.Cs(groupIdxs);
R = params.Rs(groupIdxs);
xDim = size(C{1},2);

% Compute cross-covariance matrix
Sig_12 = C{1}*C{2}';

% Compute covariant modes
[V1, P, V2] = svd(Sig_12, 'econ');
P = P(1:xDim,1:xDim);
V{1} = V1(:,1:xDim);
V{2} = V2(:,1:xDim);
