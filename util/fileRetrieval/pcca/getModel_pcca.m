function res = getModel_pcca(runIdx,xDim,varargin)
% 
% res = getModel_pcca(runIdx,xDim,...)
%
% Description: Retrieve the parameters of a pCCA model corresponding to
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
%     xDim        -- int; number of across-group latent 
%                    variables
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
%     res -- structure containing saved pCCA results (see fit_pcca.m for
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

fname = generate_filename_pcca(runIdx, xDim, ...
                               'baseDir', baseDir, 'cvf', cvf);
if ~isfile(fname)
    fprintf('ERROR: %s does not exist.  Exiting...\n', fname);
    return
else
    res = load(fname);
end
