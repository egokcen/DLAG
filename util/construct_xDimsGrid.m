function xDims_grid = construct_xDimsGrid(xDim_total)
%
% xDims_grid = construct_xDimsGrid(xDim_total)
%
% Description: Construct a grid a within- and across-group dimensionalities
%              constrained to add up to total dimensionalities for each
%              group.
%
% Arguments:
%
%     xDim_total  -- (1 x numGroups) array; total number of latents 
%                    (within and across) in each group
%
% Outputs:
%
%     xDims_grid  -- (min(xDim_total)+1 x numGroups+1) array; xDims_grid(i,:) gives
%                    the desired number of across- and within-group
%                    dimensions for model i. The first column of
%                    xDims_grid(i,:) specifies the number of across-group
%                    dimensions. The remaining columns specify the number
%                    of within-group dimensions for each group. For
%                    example,
%                        xDims_grid = [1 2 3; 4 5 6]
%                    specifies model 1 with xDim_across = 1,
%                    xDim_within = [2 3]; model 2 has xDim_across = 4, 
%                    xDim_within = [5 6].
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     13 Mar 2021 -- Initial full revision.  

numGroups = length(xDim_total);
% Largest across-group dimensionality we'll consider
xDim_across_max = min(xDim_total);
% All across-group dimensionalities we'll consider
xDims_across = 0:xDim_across_max;
% Use within-group dimensionalities such that
% xDim_across + xDim_within = xDim_total
xDims_grid = nan(length(xDims_across),numGroups+1);
for modelIdx = 1:length(xDims_across)
    xDim_across = xDims_across(modelIdx);
    xDim_within = xDim_total - xDim_across;
    xDims_grid(modelIdx,:) = [xDim_across xDim_within];
end
