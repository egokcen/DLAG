function [Xorth, Corth] = predictiveProjection_pcca(Ys, params, pDim, rGroups)
%
% [Xorth, Corth] = predictiveProjection_pcca(Ys, params, pDim, rGroups)
%
% Description: Project activity of rGroups(1) onto the subspace most 
%              predictive of rGroups(2).
%
% Arguments:
%     Ys      -- (1 x numGroups) cell array; list of data matrices 
%                {(y1Dim x N), (y2Dim x N), ...}
%     params  -- learned pCCA parameteres (structure with fields Cs, Rs, ds)
%     pDim    -- int; number of top orthonormal latent coordinates to use 
%                for projection
%     rGroups -- (1 x 2) int array; Indices of the groups we wish
%                to regress on.
%
% Outputs:
%     Xorth -- (pDim x N) array; posterior mean
%     Corth -- (yDims(rGroups(2)) x pDim) array; Top pDim predictive
%                dimensions
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Mar 2019 -- Initial full revision.

in = rGroups(1);  % Source group
out = rGroups(2); % Target group

Yin = Ys(in);
paramsReg.Cs = params.Cs(in);
paramsReg.ds = params.ds(in);
paramsReg.Rs = params.Rs(in);

% Infer latents given only the source group
[X, ~] = pcca_estep(Yin, paramsReg);

% Orthonormalize with respect to the *target* group
Ctarget = params.Cs{out};
[Xorth, Corth]   = orthogonalize(X.mean, Ctarget);

Xorth = Xorth(1:pDim,:);
Corth = Corth(:,1:pDim);
