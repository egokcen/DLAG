function [outseq, outparams] = interpolatedInference(inseq, inparams, inres, outres, varargin)
%
% [outseq, outparams] = interpolatedInference(inseq, inparams, inres, outres, ...)
%
% Description: Infer latent trajectories with resolution outres, given 
%              observations sampled with resolution (sample period or bin
%              width) inres.
%
% Arguments:
%
%     Required:
%
%     inseq     -- input data structure, whose nth entry (corresponding to
%                  the nth trial) has fields
%                      trialId        -- unique trial identifier
%                      T (1 x 1)      -- number of *input* timesteps
%                      y (yDim x Tin) -- neural data
%
%     inparams  -- Structure containing DLAG model parameters, in units
%                  normalized by input resolution, inres. 
%                  Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in units of time are given by 
%                                    'inres ./ sqrt(gamma)'                                                    
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
%                                    time-steps, spaced according to inres.
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number of
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
%     inres     -- float; input resolution (sample period or bin width), in
%                  units of time.
%     outres    -- float; output resolution, in units of time.
%
%     NOTE: inres is assumed to be a multiple of outres:
%           (inres / outres) = (integer).
%
%     Optional:
%
%     meanOnly - logical; if true, compute only posterior means, no
%                variances. Saves memory. (default: true)
%
% Outputs:
%
%     outseq    -- output data structure, whose nth entry (corresponding to
%                  the nth trial) has fields
%                  
%                  trialId        -- unique trial identifier
%                  T (1 x 1)      -- number of *output* timesteps
%                  xsm   -- ((numGroups*xDim) x Tout) array; posterior mean 
%                           at each timepoint
%                  if meanOnly == true:
%                      Vsm   -- (xDim*numGroups x xDim*numGroups x Tout) array;
%                               posterior covariance at each timepoint
%                      VsmGP_across -- (numGroups*Tout x numGroups*Tout x xDim_across)
%                                      array; posterior covariance of each 
%                                      across-group GP
%                      VsmGP_within -- (1 x numGroups) cell array;
%                                      VsmGP_within{i} -- (Tout x Tout x xDim_within(i))
%                                      array; posterior covariance of each
%                                      within-group GP for group i
%                                      VsmGP_within(i) is empty wherever 
%                                      xDim_within(i) is 0
%
%     outparams  -- Structure containing DLAG model parameters, in units
%                   normalized by output resolution, outres. 
%                   Contains the fields
% 
%                    covType -- string; type of GP covariance (e.g., 'rbf')
%                    gamma_across -- (1 x xDim_across) array; GP timescales
%                                    in units of time are given by 
%                                    'outres ./ sqrt(gamma)'                                                    
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
%                                    time-steps, spaced according to outres.
%                    xDim_across  -- int; number of across-group latent 
%                                    variables
%                    xDim_within  -- (1 x numGroups) array; number of
%                                    within-group latents in each group
%                    yDims        -- (1 x numGroups) array; 
%                                    dimensionalities of each observed group
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     26 Feb 2021 -- Initial full revision.
%     27 Feb 2021 -- Fixed indexing issue that produced future predictions.
%     13 Mar 2021 -- Added use of convertParamUnits

% Optional arguments
meanOnly = true;
assignopts(who,varargin);

% Initialize other relevant variables
yDims = inparams.yDims;
q = sum(yDims);       % shorthand for yDim
xDim_across = inparams.xDim_across;
xDim_within = inparams.xDim_within;
numGroups = length(yDims);

% Total number of within- and across-group latents, for all groups
mp = numGroups * xDim_across + sum(xDim_within);

% Precomputations for the DLAG E-step (latent inference)
precomp = makePrecomp_dlag_Estep(inseq, inparams);
CRinvC = precomp.CRinvC;

% Convert inparams to outparams with resolution outres
outparams = convertParamUnits(inparams, inres, outres);

% Initialize outseq
for n = 1:length(inseq)
    outseq(n).trialId = inseq(n).trialId;
end

% Group trials of same length together
Tall = [inseq.T];
Tu = unique(Tall);

% Overview:
% - Outer loop on each element of Tu.
% - For each element of Tu, find all trials with that length.
% - Do inference for all those trials together.

for j = 1:length(Tu)
    Tin = Tu(j); % Number of input samples
    Tout = (inres / outres) * (Tin - 1) + 1; % Number of output samples
    Tin_map = 1:(inres/outres):Tout; % Indexes of input samples at output resolution
    
    % Construct K_big based on output resolution
    [Ko] = make_K_big_dlag(outparams,Tout);
    
    % Get the subsets of K_big_out corresponding to input samples
    inBlockIdxs = get_block_idxs(repmat(mp,1,Tout));
    inIdxs = [];
    for tIdx = 1:Tin
        t = Tin_map(tIdx);
        inBlock = inBlockIdxs{t};
        inIdxs = [inIdxs inBlock(1):inBlock(2)];
    end
    Ki = Ko(inIdxs,inIdxs);
    Koi = Ko(:,inIdxs);
    
    Ki_inv = inv(Ki);
    Ko = sparse(Ko);
    
    CRinvC_big = cell(1,Tin);
    [CRinvC_big{:}] = deal(CRinvC);
    invMi = inv(Ki_inv + blkdiag(CRinvC_big{:}));
    
    if ~meanOnly
        % Posterior covariance
        KCRC = Koi * blkdiag(CRinvC_big{:});
        Sig = Ko - KCRC*Koi' + KCRC*invMi*KCRC';

        % (xDim*numGroups) X (xDim*numGroups) Posterior covariance for each timepoint
        Vsm = nan(mp,mp,Tout);
        idx = 1:mp;
        for t = 1:Tout
            cIdx = mp*(t-1)+idx;
            Vsm(:,:,t) = Sig(cIdx,cIdx);
        end

        % (numGroups*T) x (numGroups*T) Posterior covariance for each GP

        % Posterior covariance for across-group
        VsmGP_across = nan(numGroups*Tout,numGroups*Tout,xDim_across);
        idxs = extractAcrossIndices(xDim_across,xDim_within,Tout,numGroups);
        for i = 1:xDim_across
            idx = idxs{i};
            VsmGP_across(:,:,i) = Sig(idx,idx);
        end

        % Posterior covariances for within-group latents    
        VsmGP_within = cell(1,numGroups);
        idxs = extractWithinIndices(xDim_across,xDim_within,Tout,numGroups);
        for groupIdx = 1:numGroups
            % Leave VsmGP_within empty if xDim_within is 0
            if xDim_within(groupIdx) > 0
                VsmGP_within{groupIdx} = nan(Tout,Tout,xDim_within(groupIdx));
                for i = 1:xDim_within(groupIdx)
                    idx = idxs{groupIdx}{i};
                    VsmGP_within{groupIdx}(:,:,i) = Sig(idx,idx);
                end
            end
        end
    end
    
    % Process all trials with length T
    nList    = find(Tall == Tin);
    dif = precomp.Tu(j).dif;
    term1Mat = precomp.Tu(j).term1Mat;

    % Compute blkProd = CRinvC_big * invM efficiently
    blkProd = zeros(mp*Tin, mp*Tin);
    idx     = 1: mp : (mp*Tin + 1);
    for t = 1:Tin
      bIdx            = idx(t):idx(t+1)-1;
      blkProd(bIdx,:) = CRinvC * invMi(bIdx,:);
    end
    
    blkProd = Koi * (speye(mp*Tin, mp*Tin) - blkProd);   
    xsmMat  = blkProd * term1Mat; % (xDim*numGroups*Tout) x length(nList)
    
    ctr = 1;
    for n = nList
      outseq(n).T     = Tout;
      outseq(n).xsm   = reshape(xsmMat(:,ctr), mp, Tout); 
      if ~meanOnly
          outseq(n).Vsm   = Vsm;
          outseq(n).VsmGP_across = VsmGP_across;
          outseq(n).VsmGP_within = VsmGP_within;
      end
      ctr = ctr + 1;
    end
    
end
