function [result, err] = parseFilename_fa(str)
%
% [result, err] = parseFilename_fa(str)
%
% Description: Extracts method, groupIdx, xDim, and cross-validation fold 
%              from results filename, where filename has format 
%              [method]_groupIdx[groupIdx]_xDim[xDim]_cv[cvf].mat.
%
% Arguments:
%
%     str  -- string; filename to parse
%
% Outputs:
%
%     result -- structure with fields method, groupIdx, xDim, and cvf
%     err    -- boolean; indicates if input string is invalid filename
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

result = [];
err    = false;

undi = find(str == '_');
if isempty(undi)   
    err    = true;
    return
end

result.method = str(1:undi(1)-1);    
[result.groupIdx, count, errmsg] = sscanf(str(undi(1)+1:undi(2)-1), 'group%d');

[A, count, errmsg] = sscanf(str(undi(2)+1:end), 'xDim%d_cv%d.mat');

if (count < 1) || (count > 2)
    err = true;
    return
end

result.xDim = A(1);
if count == 1
    result.cvf = 0;
else
    result.cvf = A(2);
end
