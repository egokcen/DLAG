function idxs = extractAcrossIndices(xDim_across,xDim_within,T,numGroups)
%
% idxs = extractAcrossIndices(xDim_across,xDim_within,T,numGroups)
%
% Description: Extract indices into K_big or invM that correspond to the
%              across-group latents
%
% Arguments:
%
%     xDim_across  -- int; number of across-group latent 
%                     variables
%     xDim_within  -- (1 x numGroups) array; number 
%                     within-group latents in each group
%     T            -- int; number of timesteps
%     numGroups    -- int; number of groups
%
% Outputs:
%
%     idxs        -- (1 x xDim_across) cell array; 
%                    idxs(i) -- (1 x numGroups*T) array; indices 
%                               corresponding to across-group latents
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.  

    idxs = cell(1,xDim_across);
    for i = 1:xDim_across
        idxs{i} = [];
        for t = 1:T
            bigBaseIdx = (xDim_across*numGroups + sum(xDim_within))*(t-1);                
            for groupIdx = 1:numGroups
                if groupIdx > 1
                    bigIdx = bigIdx + (groupIdx-1) * xDim_across + sum(xDim_within(1:groupIdx-1));
                else
                    bigIdx = bigBaseIdx + 1;
                end  
                
                % Fill in across-group entries
                idxs{i} = [idxs{i} bigIdx+i-1];
                
            end
        end
    end

end

