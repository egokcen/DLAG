function res = learnObsParams_dlag_freq(seq, params, Xspec, varargin)
%
% learnObsParams_dlag_freq
%
% Description: Update parameters of the DLAG observation model given 
%              inferred (frequency domain) latent states.
%
% Arguments:
%
%     Required:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId         -- unique trial identifier
%                     T (1 x 1)       -- number of timesteps
%                     yfft (yDim x T) -- unitary FFT neural data
%                     xfft (xDim x T) -- unitary FFT of the latent 
%                                        posterior mean at each frequency
% 
%     params -- Structure containing DLAG model parameters. 
%               Contains the fields
%
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
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
%     Optional:
%    
%     varFloor -- float; Set private variance floor, for entries of
%                 observation noise covariance matrix. (default: 0)
%
% Outputs:
%
%     res -- Structure containing updated observation model parameters:
%            C, d, and R (see 'params' above)
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Jul 2023 -- Initial full revision.
%     13 Aug 2024 -- Added a few runtime optimizations.

% Optional arguments
varFloor = 0;
assignopts(who, varargin);

% Initialize other relevant variables
xDim_across  = params.xDim_across;
xDim_within  = params.xDim_within;
yDims        = params.yDims;   % Number of features in each group
yDim         = sum(yDims);     % Total number of features, all groups
numGroups    = length(yDims);  % Number of groups
N            = length(seq(:)); % Number of trials
Tall         = [seq.T];        % List of trial lengths
Tu           = unique(Tall);   % Group trials of same length together

% Useful for extracting correct-sized blocks from matrices later
obs_block_idxs = get_block_idxs(yDims);
lat_block_idxs = get_block_idxs([xDim_across xDim_within]);
lat_block_idxs = lat_block_idxs(2:end);
acrossIdxs = 1:xDim_across;

groupParams = partitionParams_dlag(params);

% d update
ds = cell(1,numGroups);
for groupIdx = 1:numGroups
    ds{groupIdx} = zeros(yDims(groupIdx),1); % Initialize d
    currObsGroup = obs_block_idxs{groupIdx};
    obsIdxs = currObsGroup(1):currObsGroup(2); % Current observation group
    % Collect indices for the across-group and current within-group
    % latents
    currLatGroup = lat_block_idxs{groupIdx};
    withinIdxs = currLatGroup(1):currLatGroup(2);
    latIdxs = [acrossIdxs withinIdxs];
    xDim = length(latIdxs);
    Cm = groupParams{groupIdx}.C; % Loadings for current group
    for j = 1:length(Tu)
        T = Tu(j);
        % Process all trials with length T
        nList    = find(Tall == T);
        Y = cat(3,seq(nList).yfft);
        X = cat(3,seq(nList).xfft);
        % Find the index of zero frequency
        zeroIdx = floor(T/2)+1;
        Ym = permute(Y(obsIdxs,zeroIdx,:), [1 3 2]);
        Xm = permute(X(latIdxs,zeroIdx,:), [1 3 2]);
        % Compute d update from current trial group
        ds{groupIdx} = ds{groupIdx} + sqrt(T).*sum(Ym - Cm * Xm,2);
    end
    % Final normalization factor
    ds{groupIdx} = (1/sum(Tall)).*ds{groupIdx};
end
% Collect outputs
res.d = vertcat(ds{:}); 

% C update
Cs = cell(1,numGroups);
% Intermediate terms that will be re-used for R update
YY = cell(1,numGroups);
YX = cell(1,numGroups);
XX = cell(1,numGroups);
for groupIdx = 1:numGroups
    currObsGroup = obs_block_idxs{groupIdx};
    obsIdxs = currObsGroup(1):currObsGroup(2); % Current observation group
    % Collect indices for the across-group and current within-group
    % latents
    currLatGroup = lat_block_idxs{groupIdx};
    withinIdxs = currLatGroup(1):currLatGroup(2);
    latIdxs = [acrossIdxs withinIdxs];
    xDim = length(latIdxs);
    % Initialize intermediate terms
    YY{groupIdx} = zeros(yDims(groupIdx),1);
    YX{groupIdx} = zeros(yDims(groupIdx),xDim);
    XX{groupIdx} = zeros(xDim);
    
    % Process trials of the same length
    for j = 1:length(Tu)
        T = Tu(j);
        nList = find(Tall == T);
        Y = cat(3,seq(nList).yfft);
        X = cat(3,seq(nList).xfft);

        % Handle even and odd sequence lengths
        freqs = (-floor(T/2):floor((T-1)/2))./T;
        % Find the index of zero frequency
        zeroIdx = floor(T/2)+1;

        % Zero-center observed data
        Y(obsIdxs,zeroIdx,:) = Y(obsIdxs,zeroIdx,:) - sqrt(T).*res.d(obsIdxs);
        % Variance of each neuron (used later in R update)
        YY{groupIdx} = YY{groupIdx} + sum(sum(Y(obsIdxs,:,:).*conj(Y(obsIdxs,:,:)),3),2);
        
        % Construct Qm, the time-delay operator for the current group
        if xDim_across > 0
            Qm = [exp(-1i*2*pi*params.DelayMatrix(groupIdx,:)'*freqs); ...
                  ones(xDim_within(groupIdx),T)];
        else
            Qm = ones(xDim_within(groupIdx),T);
        end

        % Inner product of Y and X
        QX = repmat(Qm,1,length(nList)) .* reshape(X(latIdxs,:,:),xDim,[]);
        YX{groupIdx} = YX{groupIdx} ...
            + reshape(Y(obsIdxs,:,:),yDims(groupIdx),[]) ...
            * QX';
        % Second moment of X
        Qm_big = repmat(permute(Qm,[1 3 2]),[1 xDim 1]);
        Qm_herm_big = permute(conj(Qm_big),[2 1 3]);
        Sm = length(nList).*(sum(Qm_big.*Xspec(j).Sx_post(latIdxs,latIdxs,:).*Qm_herm_big,3));
        XX{groupIdx} = XX{groupIdx} + Sm + QX * QX';

    end
    
    % Compute C update from intermediate terms
    Cs{groupIdx} = real(YX{groupIdx}) / real(XX{groupIdx});
    
end

% Collect outputs
res.C = blkdiag(Cs{:});

% R update
rs = cell(1,numGroups);
for groupIdx = 1:numGroups
    rs{groupIdx} = (1/sum(Tall)).*( ...
        YY{groupIdx} ...
        - 2*real(sum(YX{groupIdx}.*Cs{groupIdx},2))...
        + sum((Cs{groupIdx}*XX{groupIdx}).*Cs{groupIdx},2));
end

% Collect outputs
r = real(vertcat(rs{:}));

% Set minimum private variance
r = max(varFloor, r);
res.R = diag(r);

if any(diag(res.R) == varFloor)
    fprintf('Warning: Private variance floor used for one or more observed dimensions in DLAG.\n');
end
