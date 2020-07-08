function [estParams,seq,LL,iterTime,D,gams_across,gams_within,err_status,msg] ...
          = em_dlag(currentParams,seq,varargin)
%
% [estParams,seq,LL,iterTime,D,gams_across,gams_within] = em_dlag(currentParams,seq,...)
%
% Description: Fit DLAG model parameters using the EM algorithm.
%
% Arguments:
%
%     Required:
%
%     currentParams -- Structure containing DLAG model parameters at which EM
%                      algorithm is initialized. Contains the fields
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
%     seq      -- data structure, whose nth entry (corresponding to
%                 the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
%
%     Optional:
%
%     maxIters  -- int; number of EM iterations to run (default: 1e6)
%     tolLL     -- float; stopping criterion #1 for EM (based on LL) 
%                  (default: 1e-8)
%     tolParam  -- float; stopping criterion #2 for EM (based on delays and
%                  timescales; i.e., if across-group delays and timescales 
%                  stop changing, stop training.) (default: -Inf)
%     freqLL    -- int; data likelihood is computed every freqLL EM iterations.
%                  freqLL = 1 means that data likelihood is computed every
%                  iteration. (default: 10)
%     freqParam -- int; store intermediate values for delays and timescales
%                  and check for convergence every freqParam EM iterations 
%                  (default: 100)
%     verbose   -- logical; specifies whether to display status messages
%                  (default: true)
%     maxDelayFrac -- float in range [0,1]; Constrain estimated delays to
%                  be no more than a certain fraction of the trial length.
%                  (default: 0.5)
%     minVarFrac -- float; Set private variance floor, for entries of
%                  observation noise covariance matrix. (default: 0.01)
%     parallelize -- logical; Here, this setting just determines which 
%                    status messages will be printed. (default: false)
%     learnDelays -- logical; If set to false, then delays will remain
%                  fixed at their initial value throughout training. 
%                  (default: true)
%     learnObs  -- logical; If set to false, then observation parameters
%                  will remain fixed at their initial value throughout 
%                  training. (default: true)
%
% Outputs:
%
%     estParams -- Structure containing DLAG model parameters returned by 
%                  EM algorithm (same format as currentParams)
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
%     NOTE: The outputs below track progress throughout the fitting
%           procedure. They can be useful in debugging, if fitting produces
%           unintuitive estimates, convergence seems slow, etc.
%     LL        -- (1 x numIters) array; data log likelihood after each EM
%                  iteration
%     iterTime  -- (1 x numIters) array; computation time for each EM
%                  iteration
%     D         -- (1 x numIters) cell array; the estimated delay matrix
%                  after each EM iteration.
%     gams_across -- (1 x numIters) cell arry; estimated gamma_across after
%                    each EM iteration.
%     gams_within -- (1 x numGroups) cell arry;
%                    gams_within(i) -- (1 x numIters) cell array; estimated
%                                      gamma_within for group i after each 
%                                      EM iteration.
%     err_status -- int; 1 if data likelihood decreased during fitting. 0
%                   otherwise.
%     msg        -- string; A message indicating why fitting was stopped 
%                   (for both error and non-error cases).
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision. 
%     17 Apr 2020 -- Added 0-within-group dimension functionality
%     23 May 2020 -- Cleaned up 'verbose' functionality. Improved 
%                    exception handling if data likelihood decreases during
%                    fitting.
%     26 May 2020 -- Changed default on tolParam from 1e-4 to -Inf (i.e.,
%                    the default is to not check this criterion).
%     27 Jun 2020 -- Added 0-across-group dimension functionality
           
% Optional arguments
maxIters     = 1e6;
tolLL        = 1e-8;
tolParam     = -Inf;
freqLL       = 10;
freqParam    = 100;
verbose      = true;
maxDelayFrac = 0.5;
minVarFrac   = 0.01;
parallelize  = false;
learnDelays  = true;
learnObs     = true;
extra_opts   = assignopts(who,varargin);

% Initialize other variables
yDims        = currentParams.yDims;
yDim         = sum(yDims);
xDim_across  = currentParams.xDim_across;
xDim_within  = currentParams.xDim_within;
numGroups    = length(yDims);
LL           = [];   % Log-likelihood at each iteration
LLi          = -inf; % Initial log-likelihood
iterTime     = [];   % Time it takes to complete each iteration
deltaD_i     = inf;  % Initial change in delays between iterations
deltaGam_across_i = inf; % Initial change in across-group timescales 
                         % between iterations

% Make sure groups are valid
assert(yDim == sum(yDims));
assert(length(yDims) == length(xDim_within));

% Make sure initial delays are within specified constraints
% Convert maxDelayFrac to units of "time steps", the same units as in
% DelayMatrix
currentParams.maxDelay = round(maxDelayFrac*min([seq.T]));
maxDelay = currentParams.maxDelay;                          

% If delays are outside the range (minDelay,maxDelay), then replace them
% with a random number in [0,1].
currentParams.DelayMatrix(currentParams.DelayMatrix >= maxDelay) = rand;
currentParams.DelayMatrix(currentParams.DelayMatrix <= -maxDelay) = rand;

% Track estimated delays and timescales from iteration to iteration
D = {currentParams.DelayMatrix};
gams_across = {currentParams.gamma_across};
for groupIdx = 1:numGroups
    gams_within{groupIdx} = {currentParams.gamma_within{groupIdx}};
end

% Track error if data likelihood decreases
err_status = 0;

% Begin EM iterations
for i =  1:maxIters
    
    if verbose && ~parallelize
        fprintf('EM iteration %3d of %d ', i, maxIters);
    end
    
    % Determine when to actually compute the log-likelihood
    if (rem(i, freqLL) == 0) || (i<=2) || (i == maxIters)
        getLL = true;
    else
        getLL = false;
    end
    
    
    %% === E STEP ===
    
    if ~isnan(LLi)
        LLold   = LLi;
    end                                         
    tic; % For tracking the computation time of each iteration                              
    [seq,LLi] = exactInferenceWithLL_dlag(seq, currentParams, ...
                                          'getLL', getLL);   
                                    
    LL = [LL LLi];
    
    %% === M STEP ===
    if learnObs
        % Learn C,d,R     
        res = learnObsParams_dlag(seq, currentParams, 'minVarFrac', minVarFrac);
        currentParams.C = res.C;
        currentParams.d = res.d;
        currentParams.R = res.R;
    end
    
    % Learn GP kernel params and delays
    if currentParams.notes.learnKernelParams
        % Across- and within-group GP kernel parameters can be learned
        % independently, so it's convenient to separate them here.
        [seqAcross, seqWithin] = partitionLatents(seq, xDim_across, xDim_within);
        tempParams = currentParams;
        % Don't try to learn GP parameters when xDim_across is 0
        if xDim_across > 0
            % Across-group parameters
            tempParams.xDim = xDim_across;
            tempParams.gamma = currentParams.gamma_across;
            tempParams.eps = currentParams.eps_across;
            % learnGPparams_pluDelays performs gradient descent to learn kernel
            % parameters WITH delays
            res = learnGPparams_plusDelays(seqAcross, tempParams, ...
                                           'algorithm','em',extra_opts{:});
            switch currentParams.covType
                case 'rbf'
                    currentParams.gamma_across = res.gamma; 
                    if learnDelays
                        % Only update delays if desired. Otherwise, they will
                        % remain fixed at their initial value.
                        currentParams.DelayMatrix = res.DelayMatrix;
                    end
                    if (rem(i, freqParam) == 0) || (i == maxIters)
                        % Store current delays and timescales and compute
                        % change since last computation
                        D = [D {currentParams.DelayMatrix}];
                        gams_across = [gams_across {currentParams.gamma_across}];
                        deltaD_i = max(abs(D{end}(:) - D{end-1}(:)));
                        deltaGam_across_i = max(abs(gams_across{end} - gams_across{end-1}));
                    else
                        deltaD_i = NaN;
                        deltaGam_across_i = NaN;
                    end

            end

            % NOTE: Learning GP noise variance is currently unsupported.
            if currentParams.notes.learnGPNoise
                currentParams.eps_across = res.eps;
            end
        end
        
        % Within-group parameters: Learning these parameters is identical
        % to learning GPFA parameters, and we reuse GPFA functions here.
        for groupIdx = 1:numGroups
            % Don't try to learn GP parameters when xDim_within is 0
            if xDim_within(groupIdx) > 0
                tempParams.gamma = currentParams.gamma_within{groupIdx};
                tempParams.eps = currentParams.eps_within{groupIdx};
                % learnGPparams performs gradient descent to learn kernel
                % parameters WITHOUT delays. 
                % Reused from GPFA
                res = learnGPparams(seqWithin{groupIdx}, tempParams,...
                                    'algorithm','em',extra_opts{:});
                switch currentParams.covType
                    case 'rbf'
                        currentParams.gamma_within{groupIdx} = res.gamma;  
                        if (rem(i, freqParam) == 0) || (i == maxIters)
                            % Store current timescales
                            % We won't track the change since last computation
                            % for the within-group case, for now.
                            gams_within{groupIdx} = [gams_within{groupIdx} {res.gamma}];
                        end

                end

                % NOTE: Learning GP noise variance is currently unsupported.
                if currentParams.notes.learnGPNoise
                    currentParams.eps_within{groupIdx} = res.eps;
                end
            end
        end
    end  
    tEnd    = toc;
    iterTime = [iterTime tEnd];  % Finish tracking EM iteration time
    
    % Display the most recent likelihood that was evaluated
    if verbose && ~parallelize
        if getLL
            fprintf('       lik %f\r', LLi);
        else
            fprintf('\r');
        end
    end
    % Verify that likelihood is growing monotonically
    if i<=2
        LLbase = LLi;
    end
    if (LLi < LLold)
        err_status = 1;
        msg = sprintf('Error: Data likelihood decreased from %g to %g on iteration %d\n',...
                      LLold, LLi, i);
        break;
    elseif ((LLi-LLbase) < (1+tolLL)*(LLold-LLbase))
        % Stopping criterion #1: log-likelihood not changing
        msg = sprintf('LL has converged');
        break;
    elseif (deltaD_i < tolParam) && (deltaGam_across_i < tolParam) && (learnDelays)
        % Stopping criterion #2: Across-group delays AND timescales not changing
        % Don't use parameter changes as a convergence criterion if not
        % learning delays.
        msg = sprintf('Across-group delays and timescales have converged');
        break;
    end

end

if ~err_status
    if length(LL) < maxIters
        msg = sprintf('%s after %d EM iterations.', msg, length(LL));
    else
        msg = sprintf('Fitting stopped after maxIters (%d) was reached.', maxIters);
    end
    
end

if verbose && ~parallelize
    fprintf('%s\n', msg);
end
    
estParams = currentParams;
