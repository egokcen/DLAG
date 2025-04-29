function [seqAcross, XspecAcross, seqWithin, XspecWithin] ...
    = partitionLatents_freq(seq, Xspec, xDim_across, xDim_within)
%
% [seqAcross, XspecAcross, seqWithin, XspecWithin] ...
%     = partitionLatents_freq(seq, Xspec, xDim_across, xDim_within)
%
% Description: Helper function to extract (frequency domain) across- and 
%              within-group latents separately.
%
% Arguments:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId         -- unique trial identifier
%                     T (1 x 1)       -- number of timesteps
%                     yfft (yDim x T) -- unitary FFT neural data
%                     xfft (xDim x T) -- unitary FFT of the latent 
%                                        posterior mean at each frequency
%     Xspec     -- data structure whose jth entry, corresponding
%                  to a group of trials of the same length, has fields
%                      T       -- int; number of time steps for this
%                                 trial group
%                      Sx_post -- (xDim x xDim x T) array; posterior 
%                                 spectrum at each frequency
%                      NOTE: For DLAG, posterior covariance/spectra of X 
%                            are the same for trials of the same length.
%     xDim_across -- int; number of across-group latent 
%                    variables
%     xDim_within -- (1 x numGroups) array; number 
%                    within-group latents in each group
%
% Outputs:
%
% seqAcross -- data structure containing across-group neural trajectories:
%                  xfft -- (xDim_across x T) array; unitary FFT of the 
%                          latent posterior mean at each frequency
% XspecAcross -- data structure whose jth entry, corresponding
%                to a group of trials of the same length, has fields
%                    T       -- int; number of time steps for this
%                               trial group
%                    Sx_post -- (xDim_across x xDim_across x T) array; 
%                               posterior spectrum of across-group latents 
%                               at each frequency
% seqWithin -- (1 x numGroups) cell array; each element contains a data
%              structure, which contains within-group neural trajectories
%                  xfft  -- (xDim_within(i) x T) array; unitary FFT of the 
%                           latent posterior mean at each frequency
%              If xDim_within(i) is 0, then seqWithin(i) will be empty.
% XspecWithin -- (1 x numGroups) cell array; XspecWithin{i} is a 
%                data structure whose jth entry, corresponding to a group 
%                of trials of the same length, has fields
%                    T       -- int; number of time steps for this
%                               trial group
%                    Sx_post -- (xDim_within(i) x xDim_within(i) x T) 
%                               array; posterior spectrum of across-group 
%                               latents at each frequency
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     20 Jul 2023 -- Initial full revision.

numGroups = length(xDim_within);
N         = length(seq(:));

% Useful for extracting correct-sized blocks from matrices later
lat_block_idxs = get_block_idxs([xDim_across xDim_within]);
lat_block_idxs = lat_block_idxs(2:end);
acrossIdxs = 1:xDim_across;

% Posterior means
seqAcross = {};
seqWithin = cell(1,numGroups);
for n = 1:N
    % Across-group latent posterior means
    if xDim_across > 0
        seqAcross(n).xfft = seq(n).xfft(acrossIdxs,:);
        seqAcross(n).T = seq(n).T;
    end
    % Within-group latent posterior means
    for groupIdx = 1:numGroups
        if xDim_within(groupIdx) > 0
            currLatGroup = lat_block_idxs{groupIdx};
            withinIdxs = currLatGroup(1):currLatGroup(2);
            seqWithin{groupIdx}(n).xfft = seq(n).xfft(withinIdxs,:);
            seqWithin{groupIdx}(n).T = seq(n).T;
        end
    end
end

% Posterior covariances/spectra
XspecAcross = [];
XspecWithin = cell(1,numGroups);
for j = 1:length(Xspec)
    if xDim_across > 0
        XspecAcross(j).T = Xspec(j).T;
        XspecAcross(j).Sx_post = Xspec(j).Sx_post(acrossIdxs,acrossIdxs,:);
    end
    for groupIdx = 1:numGroups
        if xDim_within(groupIdx) > 0
            currLatGroup = lat_block_idxs{groupIdx};
            withinIdxs = currLatGroup(1):currLatGroup(2);
            XspecWithin{groupIdx}(j).T = Xspec(j).T;
            XspecWithin{groupIdx}(j).Sx_post = Xspec(j).Sx_post(withinIdxs,withinIdxs,:);
        end
    end
end