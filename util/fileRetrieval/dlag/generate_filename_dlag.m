function fname = generate_filename_dlag(runIdx, xDim_across, xDim_within, varargin)
%
% fname = generate_filename_dlag(runIdx, xDim_across, xDim_within, ...)
% 
% Description: Helper function to generate filenames where DLAG fitting 
%              results will be saved.
% 
% Arguments: 
%
%     Required:
%
%     runIdx       --  results files will be saved in ./mat_results/runXXX,
%                      where XXX is runIdx
%     xDim_across  -- int; number of across-group latent 
%                     variables
%     xDim_within  -- (1 x numGroups) array; number 
%                     within-group latents in each group
%
%     Optional:
%
%     baseDir      -- string; specifies directory in which to store
%                     mat_results. (default: '.', i.e., current directory)
%     cvf          -- int; cross-validation fold (0 means model fit to all
%                     data) (default: 0)
%     method       -- string; specifies which method was used for fitting:
%                     'dlag' or 'dlag_freq' (default: 'dlag_freq')
% 
% Outputs:
%     fname -- string; path to results file
%     
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     26 Mar 2020 -- Initial full revision.

baseDir  = '.';
cvf = 0;
method = 'dlag';
assignopts(who, varargin);
numGroups = length(xDim_within);
runDir = sprintf('%s/mat_results/run%03d', baseDir, runIdx);
fname = sprintf('%s/%s_nGroups%02d_xDimA%02d_xDimW', ...
                runDir, method, numGroups, xDim_across);
for groupIdx = 1:numGroups
    fname = sprintf('%s_%02d', fname, xDim_within(groupIdx)); 
end

if cvf > 0
    fname = sprintf('%s_cv%02d', fname, cvf);
end

fname = sprintf('%s.mat', fname);