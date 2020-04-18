function idxs = extractWithinIndices(xDim_across,xDim_within,T,numGroups)
%
% idxs = extractWithinIndices(xDim_across,xDim_within,T,numGroups)
%
% Description: Extract indices into K_big or invM that correspond to the 
%              within-group latents
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
% OUTPUTS:
%
%     idxs         -- (1 x numGroups) cell array;
%                     idxs{i} -- (1 x xDim_within(i)) cell array; 
%                                idxs{i}{j} -- (1 x T) array; indices 
%                                              corresponding to 
%                                              within-group latents
%                     If xDim_within(i) is 0, then idxs(i) will be empty.
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.  
%     17 Apr 2020 -- Added 0-within-group dimension functionality

    % Initialize output structures
    idxs = cell(1,numGroups);
    for groupIdx = 1:numGroups
        if xDim_within(groupIdx) > 0
            idxs{groupIdx} = cell(1,xDim_within(groupIdx));
        end
    end
    
    for groupIdx1 = 1:numGroups
        % Don't try to fill in idxs if xDim_within is 0
        if xDim_within(groupIdx1) > 0
            for i = 1:xDim_within(groupIdx1)
                idxs{groupIdx1}{i} = [];
                for t = 1:T
                    bigBaseIdx = (xDim_across*numGroups + sum(xDim_within))*(t-1);                
                    for groupIdx2 = 1:numGroups
                        if groupIdx2 > 1
                            bigIdx = bigIdx + (groupIdx2-1) * xDim_across + sum(xDim_within(1:groupIdx2-1));
                        else
                            bigIdx = bigBaseIdx + 1;
                        end  

                        if groupIdx2 == groupIdx1
                            % Fill in within-group entries
                            idxs{groupIdx1}{i} = [idxs{groupIdx1}{i} bigIdx+xDim_across+i-1];
                        end
                    end
                end
            end
        end
    end

end

