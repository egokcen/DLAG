function [result, err] = parseFilename_dlag(str)
%
% [result, err] = parseFilename_dlag(str)
%
% Description: Extracts method, xDim (within and across) and 
%              cross-validation fold from results filename, 
%              where filename has format 
%              [method]_xDimA[xDim_across]_xDimW[groupIdx]_[xDim_within]_cv[cvf].mat.
%
% Arguments:
%
%     str -- string; filename string to parse
%
% Outputs:
%
%     result   -- structure with fields method, xDim, and cvf
%     err      -- logical; indicates if input string is invalid filename
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     26 Mar 2020 -- Initial full revision.

  result = [];
  err    = false;
  
  undi = find(str == '_');
  if isempty(undi)   
    err    = true;
    return
  end
  
  result.method = str(1:undi(1)-1);  
  [numGroups, count, errmsg] = sscanf(str(undi(1)+1:undi(2)-1), 'nGroups%d');
  
  tmp = 'xDimA%d_xDimW';
  for groupIdx = 1:numGroups
      tmp = [tmp '_%d']; 
  end
  [A, count, errmsg] = sscanf(str(undi(2)+1:end), [tmp '_cv%d.mat']);
    
  if (count < 1) || (count > numGroups + 2)
    err = true;
    return
  end
  
  result.xDim_across = A(1);
  result.xDim_within = A(2:numGroups+1)';
  if count == numGroups+1
    result.cvf = 0;
  else
    result.cvf = A(end);
  end
