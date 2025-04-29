function precomp = makePrecomp_freq(seq, Xspec, xDim)
%
% precomp = makePrecomp_freq(seq, Xspec, xDim)
%
% Description: Precompute posterior covariances for the DLAG fitting 
%              procedure.
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
%
% Outputs
%     precomp -- Structure whose ith entry contains precomputations for 
%                the i-th state:
%                Tu   -- structure whose jth entry, corresponding to a
%                        group of trials of the same length, contains 
%                        the following:
%                        T     -- int; Length of all trials in this group
%                        numTrials -- int; Number of trials in this group
%                        XX    -- (T x 1) array; Posterior second moment
%                                 of this latent state
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     21 Jul 2023 -- Initial full revision.
%     13 Aug 2024 -- Only diagonal elements of XX are needed, so lower
%                    bound and gradients can be computed much more
%                    efficiently.

% Find unique numbers of trial lengths
Tall = [seq.T];
Tu = unique(Tall);

% Initialize output structure
for j = 1:length(Tu)
    for k = 1:xDim
        T = Tu(j);
        nList = find(Tall == T);
        precomp(k).Tu(j).T = T;
        precomp(k).Tu(j).numTrials = length(nList);
        precomp(k).Tu(j).XX  = zeros(T, 1);
    end
end

% Fill in output structure
for j = 1:length(Tu)
    T = Tu(j);
    nList = find(Tall == T);
    X = permute(cat(3,seq(nList).xfft), [2 3 1]); % (T x N x xDim)
    for k = 1:xDim
        Sx_post = reshape(Xspec(j).Sx_post(k,k,:),[],1);
        precomp(k).Tu(j).XX ...
            = precomp(k).Tu(j).numTrials.*Sx_post ...
            + sum(X(:,:,k) .* conj(X(:,:,k)), 2);
    end
end
