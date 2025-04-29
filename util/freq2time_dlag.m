function seq = freq2time_dlag(seq,params,varargin)
%
% seq = freq2time_dlag(seq,params,varargin)
%
% Description: Convert DLAG latent variables from the frequency domain
%              to the time domain.
%
% Arguments:
%
%     Required:
%
%     seq     -- data structure, whose nth entry (corresponding to
%                the nth trial) has fields
%                    trialId         -- unique trial identifier
%                    T (1 x 1)       -- number of timesteps
%                    y (yDim x T)    -- neural data
%                    xfft (xDim x T) -- latents in frequency domain
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
%     infield  -- string; name of frequency domain latents 
%                 (default: 'xfft')
%     outfield -- string; name of time domain latents 
%                 (default: 'xsm')
%
% Outputs:
%
%     seq     -- same as input structure, but with the additional field
%                    xsm (xDim x T) -- latents in time domain
%     
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     08 Jun 2023 -- Initial full revision.

infield = 'xfft';
outfield = 'xsm';
assignopts(who,varargin);

numGroups = length(params.yDims);
N = length(seq);
withinBlocks = get_block_idxs(params.xDim_within);

% Initialize output field
for n = 1:N
    seq(n).(outfield) = []; 
end

% Group trials of same length together
Tall = [seq.T];
Tu = unique(Tall);

for j = 1:length(Tu)
    T = Tu(j);
    
    % Handle even and odd sequence lengths
    freqs = (-floor(T/2):floor((T-1)/2))./T;
    
    % Process all trials with length T
    nList = find(Tall == T); 
    acrossIdxs = 1:params.xDim_across;
    for groupIdx = 1:numGroups
        Q = [];
        if params.xDim_across > 0
            % Across-group components
            Q = exp(-1i*2*pi*params.DelayMatrix(groupIdx,:)'*freqs);
        end
        % Within-group components
        Q = vertcat(Q, ones(params.xDim_within(groupIdx),T));
        withinIdxs = params.xDim_across+(withinBlocks{groupIdx}(1):withinBlocks{groupIdx}(2));
        keptIdxs = [acrossIdxs, withinIdxs];
        for n = nList
            seq(n).(outfield) = vertcat(seq(n).(outfield), ...
                real(sqrt(T).*ifft(ifftshift(Q.*seq(n).(infield)(keptIdxs,:),2),[],2)));
        end
    end
end
