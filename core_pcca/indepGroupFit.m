function estParams = indepGroupFit(Ys)
%
% estParams = indepGroupFit(Ys)
%
% Description: Fit parameters of a pCCA model assuming each group is
%              independent.
%
% Arguments: 
%     Ys -- (1 x numGroups) cell array; list of data matrices 
%           {(y1Dim x N), (y2Dim x N), ...}
%
% OUTPUTS:
%     estParams.Rs -- (1 x numGroups) cell array; Blocks of uniqueness 
%                     matrix {(y1Dim x y1Dim), (y2Dim x y2Dim), ...}
%     estParams.ds -- (1 x numGroups) cell array; List of data means 
%                     {(y1Dim x 1), (y2Dim x 1), ...}
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     23 Mar 2019 -- Initial full revision.
    
    numGroups = length(Ys);
    estParams.Rs = cell(1,numGroups);
    estParams.ds = cell(1,numGroups);
    for groupIdx = 1:numGroups
        estParams.Rs{groupIdx} = cov(Ys{groupIdx}',1); 
        estParams.ds{groupIdx} = mean(Ys{groupIdx},2);
    end