function precomp = makePrecomp_dlag_Estep(seq,params)
%
% makePrecomp_dlag_Estep
%
% Description: Precompute several values in the DLAG E-step.
%
% Arguments:
%
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
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
%                CRinv    -- (xDim x yDim) array; C' * Rinv
%                CRinvC   -- (xDim x xDim) array; C' * Rinv * C
%                Tu       -- structure whose jth entry, corresponding to a
%                            group of trials of the same length, contains 
%                            the following:
%                            dif      -- (yDim x sum(T)) array; Zero-
%                                        centered observations (y - d)
%                            term1Mat -- (xDim*numGroups*T) x length(nList)
%                                        array; An intermediate term,
%                                        CRinv * (y - d)
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.

% Initialize relevant parameters
yDims = params.yDims;
numGroups = length(yDims);

% Total number of within- and across-group latents, for all groups
mp = numGroups * params.xDim_across + sum(params.xDim_within);

% Precomputations
Rinv     = inv(params.R);    % (yDim x yDim)
Rinv     = (Rinv+Rinv') / 2; % ensure symmetry
logdet_R = logdet(params.R);

CRinv  = params.C' * Rinv;   % (xDim x yDim)
CRinvC = CRinv * params.C;   % (xDim x xDim)

% Collect outputs
precomp.Rinv = Rinv;
precomp.logdet_R = logdet_R;
precomp.CRinv = CRinv;
precomp.CRinvC = CRinvC;

% Group trials of same length together
Tall = [seq.T];
Tu = unique(Tall);

for j = 1:length(Tu)
    T = Tu(j);
    % Process all trials with length T
    nList    = find(Tall == T);
    dif      = bsxfun(@minus, [seq(nList).y], params.d); % yDim x sum(T)
    term1Mat = reshape(CRinv * dif, mp*T, []); % (xDim*numGroups*T) x length(nList)
    % Collect outputs
    precomp.Tu(j).dif = dif;
    precomp.Tu(j).term1Mat = term1Mat;    
end
