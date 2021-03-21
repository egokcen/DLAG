function outparams = getSubsetBootParams(inparams, xDims_across, xDims_within)
%
% outparams = getSubsetBootParams(inparams, xDims_across, xDims_within)
%
% Description: Return a bootstrapped GP parameter structure with only a
%              specified subset of parameters.
%
% Arguments:
%
%     Required:
%
%     inparams    -- Structure containing bootstrapped DLAG GP parameters.
%                    Contains the fields
%
%                    tau_across.upper -- (1 x xDim_across) array; 
%                                         confidence interval upper bound
%                                         on across-group GP timescales  
%                    tau_across.lower -- (1 x xDim_across) array; 
%                                         confidence interval lower bound
%                                         on across-group GP timescales 
%                    tau_within.upper -- (1 x numGroups) cell array; 
%                                        confidence interval upper bound
%                                        on GP timescales for each group
%                    tau_within.lower -- (1 x numGroups) cell array; 
%                                        confidence interval lower bound
%                                        on GP timescales for each group
%                    DelayMatrix.upper -- (numGroups x xDim_across) array;
%                                        confidence interval upper bound on 
%                                        the delay matrix, converted to 
%                                        units of time.
%                    DelayMatrix.lower -- (numGroups x xDim_across) array;
%                                        confidence interval lower bound on 
%                                        the delay matrix, converted to 
%                                        units of time.
%     xDims_across -- (1 x numDims) array; across-group state 
%                     dimensions to be retained in outparams.
%     xDims_within -- (1 x numGroups) cell array; each element is a vector 
%                     of within-group state dimensions to be retained in 
%                     outparams.
%
% Outputs:
%
%     outparams -- Structure containing subset of inparams entries.
%                  Same format as inparams, above.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     22 Sep 2020 -- Initial full revision.

numGroups = length(inparams.tau_within.upper); % Number of observation groups
xDim_across = length(inparams.tau_across.upper);
% bootstrapGPparams has an issue with collapsing the delay matrix for 1D models
if xDim_across == 1 
    inparams.DelayMatrix.upper = reshape(inparams.DelayMatrix.upper,numGroups,[]);
    inparams.DelayMatrix.lower = reshape(inparams.DelayMatrix.lower,numGroups,[]);
end

outparams.tau_across.upper = inparams.tau_across.upper(xDims_across);
outparams.tau_across.lower = inparams.tau_across.lower(xDims_across);
outparams.DelayMatrix.upper = inparams.DelayMatrix.upper(:,xDims_across);
outparams.DelayMatrix.lower = inparams.DelayMatrix.lower(:,xDims_across);

for groupIdx = 1:numGroups
    outparams.tau_within.upper{groupIdx} = inparams.tau_within.upper{groupIdx}(xDims_within{groupIdx});
    outparams.tau_within.lower{groupIdx} = inparams.tau_within.lower{groupIdx}(xDims_within{groupIdx});
end
