function precomp = makePrecomp_delays_freq(seq, Xspec, params)
%
% [precomp] = makePrecomp_delays_freq(seq, Xspec, params)
%
% Description: Precompute posterior covariances for the GP delay update
%              (frequency domain approximation).
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
%
%     Xspec     -- data structure whose jth entry, corresponding
%                  to a group of trials of the same length, has fields
%                      T       -- int; number of time steps for this
%                                 trial group
%                      Sx_post -- (xDim x xDim x T) array; posterior 
%                                 spectrum at each frequency
%                      NOTE: For DLAG, posterior covariance/spectra of X 
%                            are the same for trials of the same length.
%
%     params -- Structure containing DLAG model parameters. 
%               Contains the fields
% 
%               covType -- string; type of GP covariance (e.g., 'rbf')
%               gamma   -- (1 x xDim) array; GP timescales
%                          in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
%               eps     -- (1 x xDim) GP noise variances
%               if covType == 'sg'
%                   nu -- (1 x xDim) array; center frequencies for spectral
%                         Gaussians; convert to 1/time via 
%                         nu_across./binWidth
%               d            -- (yDim x 1) array; observation mean
%               C            -- (yDim x (numGroups*xDim)) array;
%                               mapping between low- and high-d spaces
%               R            -- (yDim x yDim) array; observation noise
%                               covariance 
%               DelayMatrix  -- (numGroups x xDim_across) array;
%                               delays from across-group latents to 
%                               observed variables. NOTE: Delays are
%                               reported as (real-valued) number of
%                               time-steps.
%               xDim    -- int; number of across-group latent variables
%               yDims        -- (1 x numGroups) array; 
%                               dimensionalities of each observed group
%
% Outputs
%     precomp -- Structure whose m-th entry contains precomputations for 
%                the m-th observation group:
%                Tu   -- structure whose jth entry, corresponding to a
%                        group of trials of the same length, contains 
%                        the following:
%                        T -- int; Length of all trials in this group
%                        numTrials -- int; Number of trials in this group
%                        XXCRinvC -- (xDim x xDim x T) array; An 
%                                    intermediate quantity to be reused.
%                        YXRinvC -- (yDims(m) x xDim x T) array; a common 
%                                   intermediate quantity to be reused  
%                NOTE: The first entry of precomp is empty because 
%                      delays to the first group are not updated.
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     24 Jul 2023 -- Initial full revision. 

% Initialize relevant variables
xDim = params.xDim;
yDims = params.yDims;
numGroups = length(yDims);
obs_block_idxs = get_block_idxs(yDims);
groupParams = partitionParams_dlag(params);

% Find unique numbers of trial lengths
Tall = [seq.T];
Tu = unique(Tall);

% Initialize output structure
for j = 1:length(Tu)
    for groupIdx = 1:numGroups
        T = Tu(j);
        nList = find(Tall == T);
        precomp(groupIdx).Tu(j).T = T;
        precomp(groupIdx).Tu(j).numTrials = length(nList);
        if groupIdx <= 1
            precomp(groupIdx).Tu(j).XXCRinvC = [];
            precomp(groupIdx).Tu(j).YXRinvC = [];
        else
            precomp(groupIdx).Tu(j).XXCRinvC = zeros(xDim,xDim,T);
            precomp(groupIdx).Tu(j).YXRinvC = zeros(yDims(groupIdx),xDim,T);
        end
    end
end

% Fill in output structure
for j = 1:length(Tu)
    T = Tu(j);
    nList = find(Tall == T);
    X = cat(3,seq(nList).xfft); % (xDim x T x N)
    Y = cat(3,seq(nList).yfft);
    % Find the index of zero frequency
    zeroIdx = floor(T/2)+1;
    % Zero-center observed data
    Y(:,zeroIdx,:) = Y(:,zeroIdx,:) - sqrt(T).*params.d;

    % XX
    XX = length(nList).*Xspec(j).Sx_post ...
        + sum(repmat(permute(X,[1 4 2 3]),1,xDim,1,1) ...
        .* repmat(permute(conj(X),[4 1 2 3]),xDim,1,1,1),4); % (xDim x xDim x T)

    for groupIdx = 2:numGroups
        currObsGroup = obs_block_idxs{groupIdx};
        obsIdxs = currObsGroup(1):currObsGroup(2); % Current observation group

        RinvC = (1./diag(groupParams{groupIdx}.R)) .* groupParams{groupIdx}.C(:,1:xDim);
        CRinvC = groupParams{groupIdx}.C(:,1:xDim)' * RinvC;

        % XXCRinvC
        precomp(groupIdx).Tu(j).XXCRinvC = permute(XX,[2 1 3]) .* repmat(CRinvC,1,1,T);

        % XY
        XY = sum(repmat(permute(X,[1 4 2 3]),1,yDims(groupIdx),1,1) ...
             .* repmat(permute(conj(Y(obsIdxs,:,:)),[4 1 2 3]),xDim,1,1,1),4); % (xDim x yDims(groupIdx) x T)

        % YXRinvC
        precomp(groupIdx).Tu(j).YXRinvC = permute(XY,[2 1 3]) .* repmat(RinvC,1,1,T);
    end
end
