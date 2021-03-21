function [P, U, V] = predictiveModes_pcca(params, varargin)
%
% [P, U, V] = predictiveModes_pcca(params, ...)
%
% Description: Compute the predictive modes between a source and target
%              group, which capture maximal predictive power from
%              source to target.
%                  P = V{1}'*Sig_11^(-0.5)*C{1}*C{2}'*V{2}
%                    = U{1}'*C{1}*C{2}'*V{2}
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
%                  analyze. Order matters: groupIdxs(1) gives the source
%                  group; groupIdxs(2) gives the target group.
%                  (default: [1 2])
%
% Outputs:
%
%     P  -- (xDim x xDim) array; diagonal matrix with the
%           cross-area predictive power along each predictive mode.
%     U  -- (sourceDim x xDim) array; predictive modes for the source 
%           group, which are uncorrelated but not necessarily orthogonal. 
%           U = Sig_11^(-.5)*V{1}
%     V  -- (1 x 2) cell array; V{i} -- (yDims(i) x xDim) array;
%           predictive modes for group i, which are orthogonal but
%           not necessarily uncorrelated.
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     20 Mar 2021 -- Initial full revision.

groupIdxs = [1 2];
assignopts(who,varargin);
numGroups = length(groupIdxs);

% Initialize outputs
V = cell(1,numGroups);

% Take only the parameters for the groups of interest
C = params.Cs(groupIdxs);
R = params.Rs(groupIdxs);
xDim = size(C{1},2);

% Compute least-squares matrix
Sig_11 = C{1}*C{1}' + R{1};
Sig_12 = C{1}*C{2}';
Corr_12 = Sig_11^(-0.5) * Sig_12;

% Compute predictive modes
[V1, P, V2] = svd(Corr_12, 'econ');
P = P(1:xDim,1:xDim);
V{1} = V1(:,1:xDim);
V{2} = V2(:,1:xDim);

U = Sig_11^(-0.5)*V{1};
