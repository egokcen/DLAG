function seq = generate_latents_dlag(params, T, N, varargin)
% 
% seq = generate_latents_dlag(params, T, N,...)
%
% Description: Generate N independent sequences of length T samples, 
%              according to a zero-mean Gaussian Process defined by the 
%              DLAG state model.
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
%     latentfield -- string; Name of data field in seq (default: 'xsm')
%
% Outputs:
%     seq -- structure whose nth entry (corresponding to the nth sequence)
%            has fields
%                trialId   -- unique trial (sequence) identifier  
%                T (1 x 1) -- number of timesteps
%                (latentfield) (p x T) -- latent sequence
%
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     09 Jan 2021 -- Initial full revision.

latentfield = 'xsm';
extraOpts = assignopts(who, varargin);
   
% If T is a scalar, then give all sequences the same length
if length(T) <= 1
    T = repmat(T,1,N);
end

for n = 1:N
    seq(n).trialId = n;
    seq(n).T = T(n);
    K_big = make_K_big_dlag(params, T(n)); % GP kernel matrix
    seq(n).(latentfield) = reshape(mvnrnd(zeros(1,size(K_big,1)), K_big),[],T(n));
end
