function res = getModel_fa(runIdx,groupIdx,xDim,varargin)
% 
% res = getModel_fa(runIdx,xDim,groupIdx...)
%
% Description: Retrieve the parameters of a FA model corresponding to
%              a particular run, observation group, and with specified 
%              latent dimensionalities.
%
% Arguments:
%
%     Required:
%
%     runIdx      -- int; results files will be loaded from 
%                    baseDir/mat_results/runXXX, where XXX is runIdx.
%                    baseDir can be specified by the user (see below)
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
%
%     res -- structure containing saved FA results (see fit_fa.m for
%            details)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

baseDir  = '.';
cvf = 0;
assignopts(who, varargin);
res = [];

fname = generate_filename_fa(runIdx, groupIdx, xDim, ...
                               'baseDir', baseDir, 'cvf', cvf);
if ~isfile(fname)
    fprintf('ERROR: %s does not exist.  Exiting...\n', fname);
    return
else
    res = load(fname);
end
