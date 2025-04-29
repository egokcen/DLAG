function [seq, LL, Xspec] = exactInferenceWithLL_freq(seq, params, varargin)
%
% [seq, LL, Xspec] = exactInferenceWithLL_freq(seq, params, ...)
%
% Description: Extract latent trajectories given DLAG model parameters,
%              using an approximate frequency domain approach.
%
% Arguments:
%
%     Required:
%
%     seq     -- data structure, whose nth entry (corresponding to
%                the nth trial) has fields
%                    trialId         -- unique trial identifier
%                    T (1 x 1)       -- number of timesteps
%                    yfft (yDim x T) -- unitary FFT of the neural data.
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
%                    xDim_within  -- (1 x numGroups) array; number of
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
%     Optional:
%
%     getLL      -- logical; specifies whether to compute data log 
%                   likelihood (default: true)
%     precompute -- logical; specifies whether to explicitly perform
%                   precomputations before performing the EM E-step.
%                   If precompute is false, then the precomp 
%                   argument needs to be specified with existing
%                   precomputations. (default: true)
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
%                   (default: {})
%
% Outputs:
%
%     seq       -- data structure with new fields (these fields are added
%                  to existing fields in the seq input argument)
%                  xfft  -- (xDim x T) array; unitary FFT of 
%                           the latent posterior mean at each frequency
%     LL        -- float; data log (pseudo) likelihood
%     Xspec     -- data structure whose jth entry, corresponding
%                  to a group of trials of the same length, has fields
%                      T       -- int; number of time steps for this
%                                 trial group
%                      Sx_post -- (xDim x xDim x T) array; posterior 
%                                 spectrum at each frequency
%                      NOTE: For DLAG, posterior covariance/spectra of X 
%                            are the same for trials of the same length.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Jul 2023 -- Initial full revision.
%     13 Aug 2024 -- Made several runtime optimizations.

% Optional arguments
getLL = true;
precompute = true;
precomp = {};
assignopts(who,varargin);

% Initialize other relevant variables
yDims = params.yDims;
q = sum(yDims);       % shorthand for yDim
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
numGroups = length(yDims);

% Total number of within- and across-group latents, for all groups
p = xDim_across + sum(xDim_within);
lat_block_idxs = get_block_idxs(repmat(p, [1 numGroups]));

% Precomputations for the DLAG E-step
if precompute
    precomp = makePrecomp_dlag_Estep_freq(seq, params);
end
% Initialize the delay operator matrix
Q = [zeros(numGroups,xDim_across) ones(numGroups,sum(xDim_within))];

% Group trials of same length together
Tall = [seq.T];
Tu = unique(Tall);
if getLL
    NT = sum(Tall);
    LL = -0.5 * q * NT * log(2*pi) - 0.5 * NT * precomp.logdet_R ...
        - 0.5 * precomp.LL_RinvY;
else
    LL = NaN;
end

% Overview:
% - Outer loop on each element of Tu.
% - For each element of Tu, find all trials with that length.
% - Do inference and LL computation for all those trials together.
for j = 1:length(Tu)
    T = Tu(j);
    nList = find(Tall == T);
    
    % Handle even and odd sequence lengths
    freqs = (-floor(T/2):floor((T-1)/2))./T;
    
    % Initialize output structures
    xfftMat = zeros(p,length(nList),T);
    Xspec(j).T = T;
    Xspec(j).Sx_post = nan(p,p,T);
    for n = nList
        seq(n).xfft = nan(p,T); 
    end    
    
    for f = 1:T
        % Construct S, the spectral density matrix
        Sx = make_S_dlag(params,freqs(f));
        Sx_inv = 1./Sx;
        logdet_Sx = sum(log(Sx));
        
        % Update Q, the time-delay operator
        Q(:,1:xDim_across) = exp(-1i*2*pi*freqs(f).*params.DelayMatrix); % numGroups x xDim
        
        % Posterior spectral density matrix
        Sx_post_inv = diag(Sx_inv);
        for groupIdx = 1:numGroups
            Sx_post_inv = Sx_post_inv ...
                + Q(groupIdx,:)' .* precomp.CRinvC{groupIdx} .* Q(groupIdx,:);  % p x p
        end
        Sx_post_inv = (Sx_post_inv + Sx_post_inv')./2; % Ensure Hermitian
        Sx_post = inv(Sx_post_inv);   % p x p
        Xspec(j).Sx_post(:,:,f) = Sx_post;
        logdet_Sx_post_inv = logdet(Sx_post_inv);

        % Posterior mean of X
        for groupIdx = 1:numGroups
            latBlock = lat_block_idxs{groupIdx};
            latIdxs = latBlock(1):latBlock(2);
            xfftMat(:,:,f) = xfftMat(:,:,f) + Q(groupIdx,:)' ...
                .* permute(precomp.Tu(j).CRinvY0(latIdxs,f,:), [1 3 2]);
        end
        xfftMat(:,:,f) = Sx_post * xfftMat(:,:,f);  % xDim x length(nList)

        if getLL
            % Compute data likelihood
            LL = LL - real(0.5 * length(nList) * (logdet_Sx + logdet_Sx_post_inv));
            % Construct full the time-delay operator matrix
            Q_full = make_delayOperator_dlag(params,freqs(f));
            LL  = LL + 0.5 * real( ...
                    sum( ...
                        ((permute(conj(precomp.Tu(j).CRinvY0(:,f,:)), [3 1 2]) * Q_full) * Sx_post) ...
                        .* (permute(precomp.Tu(j).CRinvY0(:,f,:), [3 1 2]) * conj(Q_full)), ...
                    [1 2]) ...
                );
        end
    end

    % Reshape X and collect it in seq
    ctr = 1;
    for n = nList
      seq(n).xfft(:,:) = permute(xfftMat(:,ctr,:), [1 3 2]);      
      ctr = ctr + 1;
    end

end
