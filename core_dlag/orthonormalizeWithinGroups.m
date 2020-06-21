function [seq, Corth] = orthonormalizeWithinGroups(seq, params, varargin)
%
% [seq, Corth] = orthonormalizeWithinGroups(seq, params, ...)
%
% Description: Orthonormalize the latent trajectories inferred by a DLAG
%              model according to shared variance explained within a group
%              (using within-group dimensions, across-group dimensions, or
%              both). Note that, in general, within- and across-group
%              dimensions are correlated. So, for example, orthonormalizing
%              with respect to only across-group dimensions does not
%              necessarily exclude within-group activity, and vice versa.
%
% Arguments:
%
%     Required:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId        -- unique trial identifier
%                     T (1 x 1)      -- number of timesteps
%                     xsm (xDim x T) -- latent trajectories
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
%     includeAcross -- logical; set to true to include across-group latents
%                      in the orthonormalization process
%     includeWithin -- logical; set to true to include within-group latents
%                      in the orthonormalization process
%
% Outputs:
%
%     seq     -- test data structure with new field xorth*, which contains
%                the orthonormalized trajectories for each group
%     Corth   -- (yDim x (numGroups*xDim)) array;
%                mapping between low- and high-d spaces. Each submatrix of
%                C is orthonormalized separately, for their group.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     11 Apr 2020 -- Initial full revision.
%     17 Apr 2020 -- Added 0-within-group dimension functionality
%     18 Jun 2020 -- Updated description to warn about interpreting
%                    trajectories orthonormalized with respect to only 
%                    within- or across-group dimensions.

includeAcross = true;
includeWithin = true;
assignopts(who,varargin);
assert(includeAcross || includeWithin);

yDims = params.yDims;
numGroups = length(yDims);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
xDim = sum([numGroups*xDim_across xDim_within]);
  
obs_block_idxs = get_block_idxs(yDims); % Index rows of C
lat_block_idxs = get_block_idxs(xDim_across + xDim_within); % Index cols of C

fn = 'xorth';
if includeAcross
    fn = sprintf('%s_across', fn);
end
if includeWithin
    fn = sprintf('%s_within', fn);
end

% Initialize output
for n = 1:length(seq)
    seq(n).(fn) = nan(xDim, seq(n).T);
end
Corths = cell(1,numGroups);
Xall = [seq.xsm];
% Orthonormalize latents for each group separately
for groupIdx = 1:numGroups
    % Leave Corths empty if includeAcross is false and xDim_within is 0.
    if includeAcross || xDim_within(groupIdx) > 0
        % Indices of the current group observations
        currObsBlock = obs_block_idxs{groupIdx};
        obsIdxs = currObsBlock(1):currObsBlock(2);

        % Indices of current group latents
        currLatBlock = lat_block_idxs{groupIdx};
        latIdxs = currLatBlock(1):currLatBlock(2);
        % Indices of current across-group latents
        acrossLatIdxs = latIdxs(1:xDim_across);
        % Indices of current within-group latents. Only include them if 
        % xDim_within is non-zero.
        if xDim_within(groupIdx) > 0
            withinLatIdxs = latIdxs(xDim_across+1:end);
        else
            withinLatIdxs = [];
        end
        % Determine which latents we wish to keep during orthonormalization
        keptLatIdxs = [];
        if includeAcross 
            keptLatIdxs = [keptLatIdxs acrossLatIdxs];
        end
        if includeWithin
            keptLatIdxs = [keptLatIdxs withinLatIdxs]; 
        end

        % Collect the specified within- or across-group latents and loading
        % dimensions for the current group
        X = Xall(keptLatIdxs,:);
        C = params.C(obsIdxs,keptLatIdxs);

        % Orthonormalize with respect to the specified dimensions in the 
        % current group
        [Xorth, Corth] = orthogonalize(X, C);
        seqTemp = seq;
        seqTemp = segmentByTrial(seqTemp, Xorth, fn);  

        % Add the orthonormalized latent sequence to the output structure
        for n = 1:length(seq)
            seq(n).(fn)(keptLatIdxs,:) = seqTemp(n).(fn);
        end

        Corths{groupIdx} = Corth;
    end
end

% Collect final orthonormalized loadings matrices for each group
Corth = blkdiag(Corths{:});

% Clean up output sequence
for n = 1:length(seq)
    seq(n).(fn)(any(isnan(seq(n).(fn)),2),:) = [];
end