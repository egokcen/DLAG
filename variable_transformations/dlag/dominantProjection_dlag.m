function seq = dominantProjection_dlag(seq, params, varargin)
%
% seq = dominantProjection_dlag(seq, params, ...)
%
% Description: Project latent trajectories inferred by a DLAG model onto
%              dominant modes, an orthormal basis for each group ordered
%              according to shared variance explained. The dominant modes
%              can be computed with respect to across-group latents,
%              within-group latents, or both.
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
%                      in the computation (default: true)
%     includeWithin -- logical; set to true to include within-group latents
%                      in the computation (default: true)
%
% Outputs:
%
%     seq     -- input data structure with new field
%                    xdom -- (xDim x T) array; trajectories (of all groups)
%                            projected onto dominant modes.
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
%     17 Mar 2021 -- Renamed 'dominantProjection_dlag'. Overhauled to 
%                    mirror correlativeProjection_dlag.m

includeAcross = true;
includeWithin = true;
assignopts(who,varargin);
assert(includeAcross || includeWithin);

numGroups = length(params.yDims);
xDim_across = params.xDim_across;
xDim_within = params.xDim_within;
xDim_total = xDim_across + xDim_within;
N = length(seq);
  
% Initialize output structure
fn = 'xdom';
for n = 1:N
    seq(n).(fn) = [];
end

% Compute dominant modes
[S, ~, V, ~] = dominantModes_dlag(params, ...
                                  'includeAcross', includeAcross, ...
                                  'includeWithin', includeWithin);

% Get appropriate latent trajectories
groupSeq = partitionObs(seq, xDim_total, 'datafield', 'xsm');
% Project latents for each group separately
for groupIdx = 1:numGroups
    % Extract within- and across-group latents
    [x_across, x_within] = partitionLatents_meanOnly(groupSeq{groupIdx}, ...
                            xDim_across, xDim_within(groupIdx));
                        
    % Collect the desired latent types
    X = [];
    if includeAcross
        Xacross = [x_across.xsm];
        X = [X; Xacross];
    end
    if includeWithin
       x_within = x_within{1};
       Xwithin = [x_within.xsm];
       X = [X; Xwithin];
    end
    
    % Project latents onto dominant modes
    Xdom = S{groupIdx}*V{groupIdx}'*X;
    
    % Add projected latents to output structure
    seqTemp = segmentByTrial(seq,Xdom,fn);
    for n = 1:N
        seq(n).(fn) = [seq(n).(fn); seqTemp(n).(fn)];
    end

end