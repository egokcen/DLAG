function outsig = getSubsetDelaySignificance(insig, xDims_across)
%
% outsig = getSubsetDelaySignificance(insig, xDims_across)
%
% Description: Return a delay significance structure with only entries
%              corresponding to the specified across-group latents.
%
% Arguments:
%
%     insig -- structure containing the following fields:
%              raw    -- (1 x xDim_across) array; significance of each delay
%                        evaluated on raw data, measured by decrease in 
%                        log-likelihood relative to the unaltered model.
%              upper  -- (1 x xDim_across) array; upper bound of bootstrap 
%                        confidence interval
%              lower  -- (1 x xDim_across) array; lower bound of bootstrap 
%                        confidence interval
%     xDims_across -- (1 x numDims) array; across-group state dimensions to
%                     be retained in outsig.
%
% Outputs:
%
%     outsig -- Structure containing subset of insig entries.
%               Same format as insig, above.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     14 Sep 2020 -- Initial full revision.

outsig.raw = insig.raw(xDims_across);
outsig.upper = insig.upper(xDims_across);
outsig.lower = insig.lower(xDims_across);
