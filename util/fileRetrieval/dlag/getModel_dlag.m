function res = getModel_dlag(runIdx,xDim_across,xDim_within,varargin)
% 
% res = getModel_dlag(runIdx,xDim_across,xDim_within,...)
%
% Description: Retrieve the parameters of a DLAG model corresponding to
%              a particular run, and with specified latent
%              dimensionalities.
%
% Arguments:
%
%     Required:
%
%     runIdx      -- int; results files will be loaded from 
%                    baseDir/mat_results/runXXX, where XXX is runIdx.
%                    baseDir can be specified by the user (see below)
%     xDim_across -- int; number of across-group dimensions
%     xDim_within -- (1 x numGroups) array; number of within-group 
%                    dimensions for each group
% 
%     Optional:
%
%     baseDir      -- string; specifies directory in which to store
%                     mat_results. (default: '.', i.e., current directory)
%     cvf          -- int; cross-validation fold (0 means model fit to all
%                     data) (default: 0)
% 
% Outputs:
%
%     res -- structure containing saved DLAG results (see fit_dlag.m for
%            details)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     09 Apr 2020 -- Initial full revision.

baseDir  = '.';
cvf = 0;
assignopts(who, varargin);
res = [];

fname = generate_filename_dlag(runIdx, xDim_across, xDim_within, ...
                               'baseDir', baseDir, 'cvf', cvf);
if ~isfile(fname)
    fprintf('ERROR: %s does not exist.  Exiting...\n', fname);
    return
else
    res = load(fname);
end
