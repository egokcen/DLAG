function precomp = makePrecomp_dlag_Estep_freq(seq,params)
%
% makePrecomp_dlag_Estep_freq
%
% Description: Precompute several values in the DLAG E-step frequency
%              domain implementation.
%
% Arguments:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId         -- unique trial identifier
%                     T (1 x 1)       -- number of timesteps
%                     yfft (yDim x T) -- unitary FFT of the neural data
%
%     params  -- Structure containing DLAG model parameters at which EM
%                algorithm is initialized. Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
%                    eps_across   -- (1 x xDim_across) GP noise variances
%                    gamma_within -- (1 x numGroups) cell array; 
%                                    GP timescales for each group
%                    eps_within   -- (1 x numGroups) cell array;
%                                    GP noise variances for each group
%                    d            -- (yDim x 1) array; observation mean
%                    C            -- (yDim x (numGroups*xDim)) array;
%                                    mapping between low- and high-d spaces
%                    R            -- (yDim x yDim) array; observation noise
%                                    covariance 
%                    DelayMatrix  -- (numGroups x xDim_across) array;
%                                    delays from across-group latents to 
%                                    observed variables. NOTE: Delays are
%                                    reported as (real-valued) number of
%                                    time-steps.
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number 
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
% Outputs :
%     precomp -- Structure containing the following E-step precomputations:
%                Rinv     -- (yDim x yDim) array; inverse of observation 
%                            noise covariance matrix, R
%                logdet_R -- float; log-determinant of R
%                LL_RinvY -- float; a term in the log-likelihood
%                CRinvC   -- (1 x numGroups) cell array; Each CRinvC{i}
%                            is a (xDim x xDim) array; C' * Rinv * C
%                Tu       -- structure whose jth entry, corresponding to a
%                            group of trials of the same length, contains 
%                            the following:
%                            CRinvY0  -- (numGroups*xDim x T x N)
%                                        array; An intermediate term,
%                                        CRinv * (y0)
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Jul 2023 -- Initial full revision.
%     13 Aug 2024 -- Added a few runtime optimizations, updated for
%                    compatibility with exactInferenceWithLL_freq.m.

% Initialize relevant parameters
yDims = params.yDims;
numGroups = length(yDims);
obs_block_idxs = get_block_idxs(yDims);
xDim = params.xDim_across + sum(params.xDim_within);
lat_block_idxs = get_block_idxs(repmat(xDim, [1 numGroups]));

% Precomputations
precomp.Rinv = diag(1./(diag(params.R)));    % (yDim x yDim)
precomp.logdet_R = sum(log(diag(params.R)));

% Restructure C
groupParams = partitionParams_dlag(params);
Cs = cell(1,numGroups);
for groupIdx1 = 1:numGroups
    % Get across-group matrix
    paramsAcross = getSubsetParams_dlag(groupParams{groupIdx1}, ...
                                        1:params.xDim_across, ...
                                        {[]});
    Cs{groupIdx1} = paramsAcross.C;
    
    % Collect within-group matrices
    paramsWithin = getSubsetParams_dlag(groupParams{groupIdx1}, ...
                                        [], ...
                                        {1:params.xDim_within(groupIdx1)});
    for groupIdx2 = 1:numGroups
        if groupIdx2 == groupIdx1
            Cs{groupIdx1} = horzcat(Cs{groupIdx1}, paramsWithin.C);
        else
            Cs{groupIdx1} = horzcat(Cs{groupIdx1}, ...
                zeros(yDims(groupIdx2),params.xDim_within(groupIdx2))); 
        end
    end
end

% CRinvC
CRinv = cell(1,numGroups);
precomp.CRinvC = cell(1,numGroups);
for groupIdx = 1:numGroups
    obsBlock = obs_block_idxs{groupIdx};
    obsIdxs = obsBlock(1):obsBlock(2);
    CRinv{groupIdx} = Cs{groupIdx}' ./ diag(params.R(obsIdxs,obsIdxs)).';
    precomp.CRinvC{groupIdx} = CRinv{groupIdx} * Cs{groupIdx};
end

% Group trials of same length together
Tall = [seq.T];
Tu = unique(Tall);

precomp.LL_RinvY = 0;
for j = 1:length(Tu)
    T = Tu(j);
    % Process all trials with length T
    nList    = find(Tall == T);
    % Find the index of zero frequency
    zeroIdx = floor(T/2)+1;
    % Zero-center the data
    Y0 = cat(3,seq(nList).yfft);
    Y0(:,zeroIdx,:) = Y0(:,zeroIdx,:) ...
        - repmat(sqrt(T).*params.d, [1,1,length(nList)]);
    % LL_RinvY, the log-likelihood term
    precomp.LL_RinvY = precomp.LL_RinvY ...
        + diag(precomp.Rinv).' * sum(Y0 .* conj(Y0), [2 3]);
    % CRinvY0, an intermediate term
    precomp.Tu(j).CRinvY0 = zeros(numGroups*xDim, T, length(nList));
    for groupIdx = 1:numGroups
        obsBlock = obs_block_idxs{groupIdx};
        obsIdxs = obsBlock(1):obsBlock(2);
        latBlock = lat_block_idxs{groupIdx};
        latIdxs = latBlock(1):latBlock(2);
        precomp.Tu(j).CRinvY0(latIdxs,:,:) = reshape( ...
            CRinv{groupIdx} * reshape(Y0(obsIdxs,:,:),yDims(groupIdx),[]), ...
        xDim, T, length(nList));
    end
end
