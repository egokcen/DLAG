function fname = generate_filename_fa(runIdx, groupIdx, xDim, varargin)
%
% fname = generate_filename_fa(runIdx, groupIdx, xDim, ...)
% 
% Description: Helper function to generate filenames where static
%              model fitting results will be saved.
% 
% Arguments: 
%
%     Required:
%
%     runIdx      -- results files will be saved in ./mat_results/runXXX,
%                    where XXX is runIdx
%     groupIdx    -- int; the observation group to which the FA model was
%                    fit
%     xDim        -- int; number of latent variables
%
%     Optional:
%
%     baseDir      -- string; specifies directory in which to store
%                     mat_results. (default: '.', i.e., current directory)
%     cvf          -- int; cross-validation fold (0 means model fit to all
%                     data) (default: 0)
% 
% Outputs:
%     fname -- string; path to results file
%     
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

baseDir  = '.';
cvf = 0;
assignopts(who, varargin);
method = 'fa';
runDir = sprintf('%s/mat_results/run%03d', baseDir, runIdx);
fname = sprintf('%s/%s_group%02d_xDim%02d', ...
                runDir, method, groupIdx, xDim);
if cvf > 0
    fname = sprintf('%s_cv%02d', fname, cvf);
end

fname = sprintf('%s.mat', fname);