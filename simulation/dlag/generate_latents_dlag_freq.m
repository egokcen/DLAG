function seq = generate_latents_dlag_freq(params, T, N, varargin)
% 
% seq = generate_latents_dlag_freq(params, T, N,...)
%
% Description: Generate N independent sequences of length T samples, 
%              according to a zero-mean Gaussian Process defined by the 
%              DLAG state model. Sequences are generated approximately in 
%              the frequency domain. The true covariances are thus 
%              circulant.
%
% Arguments:
%
%     Required:
%
%     params  -- Structure containing DLAG model parameters.
%                Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in units of time are given by 
%                                    'binWidth ./ sqrt(gamma)'                                                    
%                    eps_across   -- (1 x xDim_across) GP noise variances
%                    gamma_within -- (1 x numGroups) cell array; 
%                                    GP timescales for each group
%                    eps_within   -- (1 x numGroups) cell array;
%                                    GP noise variances for each group
%                    if covType == 'sg'
%                        nu_across -- (1 x xDim_across) array; center
%                                     frequencies for spectral Gaussians;
%                                     convert to 1/time via 
%                                     nu_across./binWidth 
%                        nu_within -- (1 x numGroups) cell array; 
%                                     center frequencies for each group
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
%     T         -- (1 x N) int array; T(n) gives the number of samples 
%                  (time length) for sequence n. If T is a scalar
%                  (length-1), then all sequences will be length-T.
%     N         -- int; number of sequences
%
%     Optional:
%
%     freqfield -- string; field name for frequency domain latents
%                  (default: 'xfft')
%     timefield -- string; field name for (time-delayed) time domain 
%                  latents (default: 'xsm')
%
% Outputs:
%     seq -- structure whose nth entry (corresponding to the nth sequence)
%            has fields
%                trialId   -- unique trial (sequence) identifier  
%                T (1 x 1) -- number of timesteps
%                (freqfield) -- (xDim x T) array; frequency domain latents
%                (timefield) -- (numGroups*xDim x T) array; (time-delayed)
%                               time domain latents
%
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     04 Sep 2023 -- Initial full revision.
%     13 Aug 2024 -- Updated with make_S_mdlag, which now returns only
%                    a vector of diagonal elements.
%                    Updated scaling in frequency domain so that time
%                    domain signals have unit variance.

freqfield = 'xfft';
timefield = 'xsm';
assignopts(who, varargin);

numGroups = length(params.xDim_within);
xDim = params.xDim_across + sum(params.xDim_within);
withinBlocks = get_block_idxs(params.xDim_within);
   
% If T is a scalar, then give all sequences the same length
if length(T) <= 1
    T = repmat(T,1,N);
end

% Group trials of same length together
Tu = unique(T);

% Initialize output structure
for n = 1:N
    seq(n).trialId = n;
    seq(n).T = T(n);
    seq(n).(freqfield) = nan(xDim,T(n));
end
    
% Generate all trials of the same length
for j = 1:length(Tu)
    Tj = Tu(j);
    nList = find(T == Tj);
    N = length(nList);

    % Handle even and odd sequence lengths
    freqs = (-floor(Tj/2):floor((Tj-1)/2))./Tj;

    X = zeros(xDim,Tj,N);
    for f = 1:Tj
        % Construct S, the spectral density matrix
        Sf = make_S_dlag(params,freqs(f));
        % Generate latents in the frequency domain
        % NOTE: Normally, a complex-valued normal random variable with 
        %       variance sigma would be generated according to
        %           sqrt(sigma/2) * (randn + 1i*randn)
        %       Here, we're generating a complex-valued signal in the
        %       frequency domain, and later we'll take only the real part
        %       of its time domain counterpart. Equivalently, we'll be 
        %       taking only the conjugate-symmetric portion of the
        %       frequency-domain signal. That step will throw out half of
        %       the signal's total power. We want signals in the time 
        %       domain, however, to have unit variance. Therefore, we'll 
        %       effectively double the variance in this next step, and
        %       remove the factor of sqrt(1/2). This choice will give us
        %       unit variance in the time domain.
        for k = 1:xDim
            X(k,f,:) = sqrt(Sf(k)).*(randn(1,1,N) + 1i.*randn(1,1,N));
        end
    end
    % Collect latents into output structure
    for n = 1:N
        seq(nList(n)).(freqfield) = X(:,:,n);
    end
end

% Convert latents to the time domain
seq = freq2time_dlag(seq, params, ...
                     'infield', freqfield, ...
                     'outfield', timefield);
