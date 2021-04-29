function [ccmap, ccmat] = crossCorrMap(seq, params, varargin)
%
% [ccmap, ccmat] = crossCorrMap(seq, params, ...)
%
% Description: Generate a series of cross-correlation maps between
%              across-group latents.
%
%              NOTE: All sequences (trials) must have the same length.
%
% Arguments:
%
%     Required:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId        -- unique trial identifier
%                     T (1 x 1)      -- number of timesteps
%                     (datafield) (xDim x T) -- latent trajectories
%
%     params  -- Structure containing DLAG model parameters.
%                Contains the fields
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
%     groupIdxs  -- (1 x 2) int array; Specify which pair of groups to
%                   analyze. Order matters: groupIdxs(1) gives the source
%                   group; groupIdxs(2) gives the target group.
%                   (default: [1 2])
%     corrProj   -- logical; set to true to compute cross-correlation maps
%                   through correlative modes. (default: false)
%     orth       -- logical; If corrProj is true, set to true to use
%                   orthonormal correlative modes. (default: false)
%     computeCov -- logical; set to true to compute an unnormalized 
%                   cross-covariance, instead of cross-correlation
%                   (default: false)
%
% Outputs:
%
%     ccmap -- (1 x xDim_across) cell array; ccmat{i} -- (T x 2*T-1) array;
%              cross-correlation map between a pair of groups through
%              across-group latent i ("unwrapped" version of 'ccmat').
%     ccmat -- (1 x xDim_across) cell array; ccmat{i} -- (T x T) array;
%              cross-correlation matrix between a pair of groups through
%              across-group latent i.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2021 -- Initial full revision.
%     28 Apr 2021 -- Added exception handling for xDim_across = 0 case.
%                    Added option to compute an unnormalized 'covariance'
%                    map.

groupIdxs = [1 2];
corrProj = false;
orth = false;
computeCov = false;
assignopts(who,varargin);

% Constants
N = length(seq);
T = seq(1).T;    % All trials must have the same length
numGroups = length(groupIdxs);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
xDim_total = xDim_across + xDim_within;

if xDim_across <= 0
   fprintf('crossCorrMap: xDim_across = 0. Returning empty structures ccmap, ccmat\n');
   ccmap = {};
   ccmat = {};
   return;
end

if corrProj
    % Compute cross-correlation maps through the correlative modes
    fn = 'xcorr';
    seq = correlativeProjection_dlag(seq, params, ...
                                     'groupIdxs', groupIdxs, ...
                                     'orth', orth);
    groupSeq = partitionObs(seq, repmat(xDim_across,1,numGroups),...
                            'datafield', fn);
else
    % Compute cross-correlation maps through raw latents
    fn = 'xsm';
    groupSeq = partitionObs(seq, xDim_total,'datafield', fn);
    for groupIdx = 1:length(xDim_total)
        % Take only across-group latents
        [x_across, ~] = partitionLatents_meanOnly(groupSeq{groupIdx}, xDim_across, xDim_within(groupIdx));
        for n = 1:N
            groupSeq{groupIdx}(n).xsm = x_across(n).xsm;
        end
    end
    groupSeq = groupSeq(groupIdxs); % Take only the pair of groups we want
end

% Actual cross-correlation map computation
% Stack all time points into one long column vector
X = cell(1,xDim_across);
for j = 1:xDim_across
    for n = 1:N
        xGrouped = [];
        for groupIdx = 1:numGroups
            x = groupSeq{groupIdx}(n).(fn)(j,1:T)';
            xGrouped = [xGrouped; x];
        end
        X{j} = [X{j} xGrouped];
    end
end

% Cross-correlation matrix
ccmat = cell(1,xDim_across);
for j = 1:xDim_across
    if computeCov
        ccmat{j} = cov(X{j}');        % Auto- and cross-covariance
    else
        ccmat{j} = corr(X{j}');       % Auto- and cross-correlation
    end
    ccmat{j} = ccmat{j}(1:T,T+1:end); % Take only cross-interactions
end

% Unwrap the cross-correlation matrix to get a cross-correlation map 
ccmap = cell(1,xDim_across);
for j = 1:xDim_across
    ccmap{j} = nan(T,2*T-1);
    ccmat_ud = flipud(ccmat{j});
    for t = 0:T-1
        ccmap{j}(1+t,1+t:t+T) = ccmat_ud(1+t,:);
    end
end
