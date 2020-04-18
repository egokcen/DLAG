function [result, err] = parseFilename_pcca(str)
%
% [result, err] = parseFilename_pcca(str)
%
% Description: Extracts method, xDim, and cross-validation fold from 
%              results filename, where filename has format 
%              [method]_xDim[xDim]_cv[cvf].mat.
%
% Arguments:
%
%     str  -- string; filename to parse
%
% Outputs:
%
%     result -- structure with fields method, xDim, and cvf
%     err    -- boolean; indicates if input string is invalid filename
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     09 Apr 2020 -- Initial full revision.

result = [];
err    = false;

undi = find(str == '_');
if isempty(undi)   
err    = true;
return
end

result.method = str(1:undi(1)-1);    

[A, count, errmsg] = sscanf(str(undi(1)+1:end), 'xDim%d_cv%d.mat');

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
