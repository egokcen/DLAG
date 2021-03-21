function fname = generate_inference_fname_dlag(runIdx, prefix, xDim_across, xDim_within, varargin)
%
% fname = generate_inference_fname_dlag(runIdx, prefix, xDim_across, xDim_within, ...)
% 
% Description: Helper function to generate filenames where DLAG inference 
%              results will be saved.
% 
% Arguments: 
%
%     Required:
%
%     runIdx       -- int; results files will be saved in ./mat_results/runXXX,
%                     where XXX is runIdx
%     prefix       -- string; additional prefix to the filename
%     xDim_across  -- int; number of across-group latent 
%                     variables
%     xDim_within  -- (1 x numGroups) array; number 
%                     within-group latents in each group
%
%     Optional:
%
%     baseDir      -- string; specifies directory in which to store
%                     mat_results. (default: '.', i.e., current directory)
% 
% Outputs:
%     fname -- string; path to results file
%     
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     20 Mar 2021 -- Initial full revision.

baseDir  = '.';
assignopts(who, varargin);

numGroups = length(xDim_within);

fname = sprintf('%s/mat_results/run%03d', baseDir, runIdx);
fname = sprintf('%s/%s_xDimA%02d_xDimW', fname, prefix, xDim_across);
for groupIdx = 1:numGroups
    fname = sprintf('%s_%02d', fname, xDim_within(groupIdx)); 
end
fname = sprintf('%s.mat', fname);