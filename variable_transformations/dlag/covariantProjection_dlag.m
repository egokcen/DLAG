function seq = covariantProjection_dlag(seq, params, varargin)
%
% seq = covariantProjection_dlag(seq, params, ...)
%
% Description: Project the latent trajectories inferred by a DLAG
%              model onto covariant modes, which capture the maximal 
%              (0-lag) covariance between a pair of groups.
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
%     groupIdxs -- (1 x 2) int array; Specify which pair of groups to
%                  analyze. Order does not change the computation, but it
%                  does change the order of groups in the outputs. 
%                  (default: [1 2])
%     zerolag   -- logical; set true to compute zero-lag modes, false
%                  to compute modes that factor in delays (default: true)
%
% Outputs:
%
%     seq     -- input data structure with new field
%                xcov  -- (2*xDim_across x T) array; trajectories
%                         projected onto covariant modes. The first 
%                         xDim_across latents correspond to group 
%                         groupIdxs(1). The next xDim_across latents 
%                         correspond to group groupIdxs(2).
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     26 Feb 2022 -- Initial full revision.
%     06 Apr 2022 -- Added zerolag option.

groupIdxs = [1 2];
zerolag = true;
assignopts(who,varargin);

% Constants
N = length(seq);
numGroups = length(groupIdxs);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
xDim_total = xDim_across + xDim_within;

% Initialize output structure
for n = 1:N
    seq(n).xcov = []; 
end

if xDim_across <= 0
   fprintf('covariantProjection_dlag: xDim_across = 0. Returning empty field seq.xcov\n');
   return;
end

% Compute covariant modes
[~, V, ~] = covariantModes_dlag(params, 'groupIdxs', groupIdxs, 'zerolag', zerolag);

% Get across-group loading matrices
groupParams = partitionParams_dlag(params);
groupParams = groupParams(groupIdxs);
Ca = cell(1,numGroups);
for groupIdx = 1:numGroups
    acrossParams = getSubsetParams_dlag(groupParams{groupIdx}, 1:xDim_across, cell(1,numGroups));
    Ca{groupIdx} = acrossParams.C; 
end

% Take only trajectories from the groups in groupIdxs
groupSeq = partitionObs(seq, xDim_total, 'datafield', 'xsm');
groupSeq = groupSeq(groupIdxs);
xDim_within = xDim_within(groupIdxs);
for groupIdx = 1:numGroups
    % Collect across-group latents
    [x_across, ~] = partitionLatents_meanOnly(groupSeq{groupIdx}, ...
                        xDim_across, xDim_within(groupIdx));
    Xacross = [x_across.xsm];
    % Project onto orthogonal basis
    Xcov = V{groupIdx}' * Ca{groupIdx} * Xacross;
    % Add projected latents to output structure
    seqTemp = segmentByTrial(seq,Xcov,'xcov');
    for n = 1:N
        seq(n).xcov = [seq(n).xcov; seqTemp(n).xcov];
    end
end
