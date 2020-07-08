function [seqAcross, seqWithin] = partitionLatents_meanOnly(seq, xDim_across, xDim_within, varargin)
%
% [seqAcross, seqWithin] = partitionLatents_meanOnly(seq, xDim_across, xDim_within, ...)
%
% Description: Helper function to extract across- and within-group latents
%              separately. Compared to partitionLatents.m,
%              partitionLatents_meanOnly.m extracts only posterior means.
%
% Arguments:
%
%     Required:
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
%     Optional
%
%     xspec       -- string; partition latents with field name given by
%                    xspec (default: 'xsm')
%
% Outputs:
%
% seqAcross -- data structure containing across-group neural trajectories:
%               (xspec)  -- ((numGroups*xDim_across) x T) array;
%                           posterior mean at each timepoint
%              If xDim_across is 0, then seqAcross will be empty.
% seqWithin -- (1 x numGroups) cell array; each element contains a data
%              structure, which contains within-group neural trajectories
%               (xspec)  -- (xDim_within(i) x T) array; posterior mean at
%                           each timepoint
%              If xDim_within(i) is 0, then seqWithin(i) will be empty.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.   
%     17 Apr 2020 -- Added 0-within-group dimension functionality
%     27 Jun 2020 -- Patched potential indexing error (acrossIdxs_start)
%                    for numbers of groups greater than 2. Added
%                    0-across-group dimension functionality.


    xspec     = 'xsm';
    extraOpts = assignopts(who, varargin);
    numGroups = length(xDim_within);
    N         = length(seq(:));
    
    seqAcross = {};
    seqWithin = cell(1,numGroups);
    for n = 1:N
        % Initialize output structures
        if xDim_across > 0
            seqAcross(n).(xspec) = [];
            seqAcross(n).T = seq(n).T;
        end
        for groupIdx = 1:numGroups
            if groupIdx > 1
                acrossIdxs_start = 1+(groupIdx-1)*xDim_across+sum(xDim_within(1:groupIdx-1));
            else
                acrossIdxs_start = 1;
            end
            acrossIdxs_end = acrossIdxs_start + xDim_across - 1;
            if xDim_across > 0
                seqAcross(n).(xspec) = [seqAcross(n).(xspec); seq(n).(xspec)(acrossIdxs_start:acrossIdxs_end,:)];
            end
            if xDim_within(groupIdx) > 0
                withinIdxs_start = acrossIdxs_end + 1;
                withinIdxs_end = withinIdxs_start + xDim_within(groupIdx) - 1;
                seqWithin{groupIdx}(n).(xspec) = seq(n).(xspec)(withinIdxs_start:withinIdxs_end,:);
                seqWithin{groupIdx}(n).T = seq(n).T;
            end
        end
    end

end
