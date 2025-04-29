function [seq, LL] = exactInferenceWithLL_dlag(seq, params, varargin)
%
% [seq, LL] = exactInferenceWithLL_dlag(seq, params, ...)
%
% Description: Extract latent trajectories given DLAG model parameters.
%
% Arguments:
%
%     Required:
%
%     seq     -- data structure, whose nth entry (corresponding to
%                the nth trial) has fields
%                    trialId      -- unique trial identifier
%                    T (1 x 1)    -- number of timesteps
%                    y (yDim x T) -- neural data
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
%     precomp    -- Structure containing the following E-step precomputations:
%                   Rinv     -- (yDim x yDim) array; inverse of observation 
%                               noise covariance matrix, R
%                   logdet_R -- float; log-determinant of R
%                   CRinv    -- (xDim x yDim) array; C' * Rinv
%                   CRinvC   -- (xDim x xDim) array; C' * Rinv * C
%                   Tu       -- structure whose jth entry, corresponding to a
%                               group of trials of the same length, contains 
%                               the following:
%                               dif      -- (yDim x sum(T)) array; Zero-
%                                           centered observations (y - d)
%                               term1Mat -- (xDim*numGroups*T) x length(nList)
%                                           array; An intermediate term,
%                                           CRinv * (y - d)
%                   (default: {})
%
% Outputs:
%
%     seq       -- data structure with new fields (these fields are added
%                  to existing fields in the seq input argument)
%                  xsm   -- ((numGroups*xDim) x T) array; posterior mean 
%                           at each timepoint
%                  Vsm   -- (xDim*numGroups x xDim*numGroups x T) array;
%                           posterior covariance at each timepoint
%                  VsmGP_across -- (numGroups*T x numGroups*T x xDim_across)
%                                  array; posterior covariance of each 
%                                  across-group GP
%                  VsmGP_within -- (1 x numGroups) cell array;
%                                  VsmGP_within{i} -- (T x T x xDim_within(i))
%                                  array; posterior covariance of each
%                                  within-group GP for group i
%                                  VsmGP_within(i) is empty wherever 
%                                  xDim_within(i) is 0
%     LL        -- float; data log likelihood
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision.   
%     17 Apr 2020 -- Added 0-within-group dimension functionality

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
mp = numGroups * xDim_across + sum(xDim_within);

% Precomputations for the DLAG E-step
if precompute
    precomp = makePrecomp_dlag_Estep(seq, params);
end
Rinv = precomp.Rinv;
logdet_R = precomp.logdet_R;
CRinv = precomp.CRinv;
CRinvC = precomp.CRinvC;

% Group trials of same length together
Tall = [seq.T];
Tu = unique(Tall);
LL = 0;

% Overview:
% - Outer loop on each element of Tu.
% - For each element of Tu, find all trials with that length.
% - Do inference and LL computation for all those trials together.
for j = 1:length(Tu)
    T = Tu(j);
    
    % Construct K_big
    [K_big] = make_K_big_dlag(params,T);
    K_big_inv = inv(K_big);
    logdet_K_big = logdet(K_big);
    K_big = sparse(K_big);
    
    CRinvC_big = cell(1,T);
    [CRinvC_big{:}] = deal(CRinvC);
    invM = inv(K_big_inv + blkdiag(CRinvC_big{:}));
    logdet_M = logdet(K_big_inv + blkdiag(CRinvC_big{:}));
    
    % (xDim*numGroups) X (xDim*numGroups) Posterior covariance for each timepoint
    Vsm = nan(mp,mp,T);
    idx = 1:mp;
    for t = 1:T
        cIdx = mp*(t-1)+idx; %idx(t):idx(t+1)-1;
        Vsm(:,:,t) = invM(cIdx,cIdx);
    end
    
    % (numGroups*T) x (numGroups*T) Posterior covariance for each GP
    
    % Posterior covariance for across-group
    VsmGP_across = nan(numGroups*T,numGroups*T,xDim_across);
    idxs = extractAcrossIndices(xDim_across,xDim_within,T,numGroups);
    for i = 1:xDim_across
        idx = idxs{i};
        VsmGP_across(:,:,i) = invM(idx,idx);
    end
    
    % Posterior covariances for within-group latents    
    VsmGP_within = cell(1,numGroups);
    idxs = extractWithinIndices(xDim_across,xDim_within,T,numGroups);
    for groupIdx = 1:numGroups
        % Leave VsmGP_within empty if xDim_within is 0
        if xDim_within(groupIdx) > 0
            VsmGP_within{groupIdx} = nan(T,T,xDim_within(groupIdx));
            for i = 1:xDim_within(groupIdx)
                idx = idxs{groupIdx}{i};
                VsmGP_within{groupIdx}(:,:,i) = invM(idx,idx);
            end
        end
    end
    
    % Process all trials with length T
    nList    = find(Tall == T);
    dif = precomp.Tu(j).dif;
    term1Mat = precomp.Tu(j).term1Mat;

    % Compute blkProd = CRinvC_big * invM efficiently
    blkProd = zeros(mp*T, mp*T);
    idx     = 1: mp : (mp*T + 1);
    for t = 1:T
      bIdx            = idx(t):idx(t+1)-1;
      blkProd(bIdx,:) = CRinvC * invM(bIdx,:);
    end
    
    blkProd = K_big * (speye(mp*T, mp*T) - blkProd);   
    xsmMat  = blkProd * term1Mat; % (xDim*numGroups*T) x length(nList)
    
    ctr = 1;
    for n = nList
      seq(n).xsm   = reshape(xsmMat(:,ctr), mp, T);      
      seq(n).Vsm   = Vsm;
      seq(n).VsmGP_across = VsmGP_across;
      seq(n).VsmGP_within = VsmGP_within;

      ctr = ctr + 1;
    end
    
    if getLL
      % Compute data likelihood
      val = -T * logdet_R - logdet_K_big - logdet_M -...
            q * T * log(2*pi);      
      LL  = LL + length(nList) * val - sum(sum((Rinv * dif) .* dif)) +...
            sum(sum((term1Mat' * invM) .* term1Mat')); 
    end
    
end

if getLL
    LL = LL / 2;
else
    LL = NaN;
end
