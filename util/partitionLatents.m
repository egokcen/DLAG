function [seqAcross, seqWithin] = partitionLatents(seq, xDim_across, xDim_within)
%
% [seqAcross, seqWithin] = partitionLatents(seq, xDim_across, xDim_within)
%
% Description: Helper function to extract across- and within-group latents
%              separately.
%
% Arguments:
%
%     seq         -- data structure, whose nth entry (corresponding to
%                    the nth trial) has fields
%                        trialId      -- unique trial identifier
%                        T (1 x 1)    -- number of timesteps
%                        y (yDim x T) -- neural data
%     xDim_across -- int; number of across-group latent 
%                    variables
%     xDim_within -- (1 x numGroups) array; number 
%                    within-group latents in each group
%
% Outputs:
%
% seqAcross -- data structure containing across-group neural trajectories:
%                  xsm   -- ((numGroups*xDim_across) x T) array;
%                           posterior mean at each timepoint
%                  VsmGP -- (numGroups*T x numGroups*T x xDim_across) array;
%                           posterior covariance of each GP
% seqWithin -- (1 x numGroups) cell array; each element contains a data
%              structure, which contains within-group neural trajectories
%                  xsm   -- (xDim_within(i) x T) array; posterior mean at
%                           each timepoint
%                  VsmGP -- (T x T x xDim_within(i)) array; posterior
%                           covariance of each GP
%              If xDim_within(i) is 0, then seqWithin(i) will be empty.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.   
%     17 Apr 2020 -- Added 0-within-group dimension functionality


    numGroups = length(xDim_within);
    N         = length(seq(:));
    
    seqWithin = cell(1,numGroups);
    for n = 1:N
        % Initialize output structures
        seqAcross(n).xsm = [];
        seqAcross(n).VsmGP = seq(n).VsmGP_across;
        seqAcross(n).T = seq(n).T;
        for groupIdx = 1:numGroups
            if groupIdx > 1
                acrossIdxs_start = 1+(groupIdx-1)*xDim_across+sum(xDim_within(1:numGroups-1));
            else
                acrossIdxs_start = 1;
            end
            acrossIdxs_end = acrossIdxs_start + xDim_across - 1;
            seqAcross(n).xsm = [seqAcross(n).xsm; seq(n).xsm(acrossIdxs_start:acrossIdxs_end,:)];
            if xDim_within(groupIdx) > 0
                withinIdxs_start = acrossIdxs_end + 1;
                withinIdxs_end = withinIdxs_start + xDim_within(groupIdx) - 1;
                seqWithin{groupIdx}(n).xsm = seq(n).xsm(withinIdxs_start:withinIdxs_end,:);
                seqWithin{groupIdx}(n).VsmGP = seq(n).VsmGP_within{groupIdx};
                seqWithin{groupIdx}(n).T = seq(n).T;
            end
        end
    end

end
