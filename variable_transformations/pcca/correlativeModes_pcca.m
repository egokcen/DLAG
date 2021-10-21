function [P, U, V] = correlativeModes_pcca(params, varargin)
%
% [P, U, V] = correlativeModes_pcca(params, ...)
%
% Description: Compute the correlative modes between a pair of groups, 
%              which capture the maximal 0-lag correlation across groups.
%                  P = V{1}'*Sig_11^(-0.5)*C{1}*C{2}'*Sig_22^(-0.5)*V{2}
%                    = U{1}'*C{1}*C{2}'*U{2}
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
%     U  -- (1 x 2) cell array; U{i} -- (yDims(groupIdxs(i)) x xDim) array;
%           correlative modes for group groupIdxs(i), which are 
%           uncorrelated but not necessarily orthogonal. 
%           U{i} = Sig_ii^(-.5)*V{i}
%     V  -- (1 x 2) cell array; V{i} -- (yDims(groupIdxs(i)) x xDim) array;
%           correlative modes for group groupIdxs(i), which are orthogonal
%           but not necessarily uncorrelated.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     20 Mar 2021 -- Initial full revision.
%     20 Oct 2021 -- Updated documentation to clarify group order in U, V.

groupIdxs = [1 2];
assignopts(who,varargin);
numGroups = length(groupIdxs);

% Initialize outputs
U = cell(1,numGroups);
V = cell(1,numGroups);

% Take only the parameters for the groups of interest
C = params.Cs(groupIdxs);
R = params.Rs(groupIdxs);
xDim = size(C{1},2);

% Compute cross-correlation matrix
Sig_11 = C{1}*C{1}' + R{1};
Sig_22 = C{2}*C{2}' + R{2};
Sig_12 = C{1}*C{2}';
Corr_12 = Sig_11^(-0.5) * Sig_12 * Sig_22^(-0.5);

% Compute correlative modes
[V1, P, V2] = svd(Corr_12, 'econ');
P = P(1:xDim,1:xDim);
V{1} = V1(:,1:xDim);
V{2} = V2(:,1:xDim);

U{1} = Sig_11^(-0.5)*V{1};
U{2} = Sig_22^(-0.5)*V{2};
