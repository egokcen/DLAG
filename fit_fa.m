function result = fit_fa(runIdx, dat, varargin)
%
% result = fit_fa(runIdx, dat,...)
%
% Description: This function does most of the front-end work for fitting
%              factor analysis models. It prepares directories for saved 
%              results, sets up cross-validation folds (along with 
%              parallelization), and then calls the functions that actually
%              estimate model parameters.
%
% Arguments: 
%
%     Required:
%
%     runIdx    --  int; results files will be saved in 
%                   baseDir/mat_results/runXXX, where XXX is runIdx.
%                   baseDir can be specified by the user (see below)
%     dat       --  (1 x N) structure whose nth entry (corresponding to the
%                   nth experimental trial) has fields
%                       trialId      -- unique trial identifier  
%                       T (1 x 1)    -- number of timesteps
%                       y (yDim x T) -- continuous valued data 
%                                       (e.g., binned spike counts)
%                   Features (neurons) are assumed to be ordered according
%                   to which group (area) they belong to. For example, if
%                   data has features belonging to two groups, then each 
%                   observation should be organized as follows: 
%                        [      |
%                         Group 1 Features
%                               |
%                               |
%                         Group 2 Features
%                               |         ]
%                   NOTE: If you have spiking neural data, and need to
%                         convert it to binned spike counts, see getSeq.m.
%
%     Optional: (Many of these are technically defined as optional, but
%               should really be specified).
%
%     baseDir     -- string; specifies directory in which to store
%                    mat_results. (default: '.', i.e., current directory)
%     method      -- string; method to be used (right now, just 'fa') 
%                    (default: 'fa')
%     binWidth    -- float; bin width or sample period, in units of time
%                    (e.g., ms) (default: 20)
%     numFolds    -- int; number of cross-validation folds (default: 0).
%                    0 indicates no cross-validation, i.e. train on all trials.
%     yDims       -- (1 x numGroups) array; Specify the number of features 
%                    (neurons) in each group (area). Elements in yDims
%                    should match the format of dat.
%     xDims       -- (1 x numGroups) cell array; 
%                    xDims{i} -- (1 x numDims) array; latent 
%                      dimensionalities to be modeled for group i. For 
%                      example, to fit 3 separate models, each with 1, 2, 
%                      and 3 dimensions, input 1:3 or [1 2 3]. 
%                      (default: {})
%     parallelize -- logical; Set to true to use Matlab's parfor construct
%                    to parallelize each fold and latent dimensionality 
%                    using multiple cores. (default: false)
%     numWorkers  -- int; Number of cores to use, if using the parallelize
%                    option. (default: 4)
%     randomSeed  -- int (or struct); Specify a seed (or full settings) for
%                    the random number generator, to aid reproducibility.
%                    For example, specifying the same seed for different 
%                    runs will generate the same cross-validation folds. 
%                    If empty ([]), seed the random number generator with 
%                    the current time (rng('shuffle')). (default: 0)
%
% Outputs:
%
%     Results are saved to a file. These variables should help with 
%     reproducibility, and include:
%               
%               fname       -- string; path to file where variables are 
%                              stored
%               method      -- string; method that was used ('fa')
%               rngSettings -- structure with the random number generator 
%                              settings used during run time. Includes 
%                              fields 'Type', 'Seed', and 'State'.
%               groupIdx    -- int; indicates the group to which this FA
%                              model was fit.
%               xDim        -- int; number of latent dimensions
%               yDim        -- observation dimensionality of the group to
%                              which this FA model was fit.
%               binWidth    -- float; bin width or sample period, in units 
%                              of time (e.g., ms)
%               parallelize -- logical; indicates whether or not 
%                              parallelization was used.
%               saveData    -- logical; indicates whether or not train and 
%                              test data is saved in this file.
%               startParams -- structure containing parameter values
%                              at which FA EM procedure was initialized.
%               estParams   -- structure containing final static model 
%                              parameter estimates, after model fitting was
%                              completed.
%               LL          -- (1 x numIters) array; data log likelihood  
%                              after each FA EM iteration
%               LLtrain     -- float; final data log likelihood evaluated 
%                              on train data
%               cvf         -- int; indicates which cross-validation fold 
%                              this model was fit to. (0 if fit to all 
%                              data)
%               seqTrain    -- training data structure of the same format 
%                              as dat. Also includes latent trajectories 
%                              inferred using estParams.
%                              NOTE: Only included if saveData is true. 
%               extraOpts   -- cell array; Contains additional optional 
%                              arguments specified by the user at runtime.
%                              These arguments may have overridden defaults.
%
%     Below are additional variables saved to files corresponding to
%     cross-validation folds (these variables will NOT appear in result
%     above):
%     
%     LLtest  -- float; final data log likelihood evaluated 
%                on test data
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 Aug 2020 -- Initial full revision.

% Specify defaults for optional arguments
baseDir              = '.';
method               = 'fa';
binWidth             = 20;
numFolds             = 0;
yDims                = [];
xDims                = {};
parallelize          = false;
numWorkers           = 4;
randomSeed           = 0;
extraOpts            = assignopts(who, varargin);
numGroups            = length(yDims); % Number of groups (areas)

fprintf('\n---------------------------------------\n');
runDir = sprintf('%s/mat_results', baseDir);
if ~isfolder(runDir)
    mkdir(runDir);
end

% Make a directory for this runIdx if it doesn't already exist
runDir = sprintf('%s/run%03d', runDir, runIdx);
if isfolder(runDir)
    fprintf('Using existing directory %s...\n', runDir);
else
    fprintf('Making directory %s...\n', runDir);
    mkdir(runDir);
end

% Verify that yDims was specified correctly
if sum(yDims) ~= size(dat(1).y,1)
    fprintf('Error: Entries in yDims do not sum to yDim.  Exiting.\n');
    return;
end

N    = length(dat);      % Number of trials
cvf_list = 0:numFolds;   % Cross-validation folds, including training on all data

% Seed the random number generator, for reproducibility
if ~isempty(randomSeed)
    rng(randomSeed);
else
    % If no seed given, then use the current time to seed the generator
    rng('shuffle');
end
% Save the random number generator settings, for reproducibility
rngSettings = rng;

% Set cross-validation folds (crossvalind randomly generates indices)
if numFolds > 0
    val_indices = crossvalind('Kfold', N, numFolds);
end

i = 1;
for cvf = cvf_list
    % Set cross-validation folds
    val = false(1, N);
    if cvf > 0
        val = (val_indices == cvf);
    end
    train = ~val;

    seqTrain      = dat(train);
    seqTest       = dat(val);
    
    % Split observations according to their groups.
    groupSeqTrain = partitionObs(seqTrain, yDims);
    groupSeqTest = partitionObs(seqTest, yDims);
    
    % Setup FA models for each group separately
    for groupIdx = 1:numGroups
        curr_xDims = xDims{groupIdx}; % List of candidate dimensionalities
        
        % Grab only observations corresponding to the current group.
        currTrain = groupSeqTrain{groupIdx};
        currTest = groupSeqTest{groupIdx};
        
        for xDim = curr_xDims
            % Check if training data covariance is full rank
            yAll = [currTrain.y];
            yDim  = size(yAll, 1);
            if rank(cov(yAll')) < yDim
                fprintf('ERROR: Observation covariance matrix for group %d is rank deficient.\n', groupIdx);
                fprintf('Possible causes: repeated units, not enough observations.\n');
                fprintf('Exiting...\n');
                return
            end

            cvf_params(i).seqTrain = currTrain;
            cvf_params(i).seqTest = currTest;

            % Specify filename where results will be saved
            fname = sprintf('%s/%s_group%02d_xDim%02d', runDir, method, groupIdx, xDim);
            if cvf > 0
                fname = sprintf('%s_cv%02d', fname, cvf);
            end
            cvf_params(i).fname = fname;
            cvf_params(i).groupIdx = groupIdx;
            cvf_params(i).xDim = xDim;
            cvf_params(i).cvf = cvf;
            i = i+1;
        end
    end
end
    
% Set up parallelization, if desired
if parallelize
    StartParPool(numWorkers);  % Helper function to set up parfor construct
    parfor i = 1:length(cvf_params)
        % Print minimal info about the model to be fitted.
        fprintf('Group: %d of %d, xDim: %d, CV Fold: %d of %d\n', ...
             cvf_params(i).groupIdx, numGroups, cvf_params(i).xDim, cvf_params(i).cvf, numFolds);
        
        % Call the FA engine
        call_fa_engine(cvf_params(i).fname, cvf_params(i).seqTrain, ...
            cvf_params(i).seqTest, 'method', method, ...
            'xDim', cvf_params(i).xDim, 'cvf', cvf_params(i).cvf, ...
            'groupIdx', cvf_params(i).groupIdx, 'parallelize', parallelize, ...
            'binWidth', binWidth, 'rngSettings', rngSettings, extraOpts{:});  
    end
else % Otherwise, continue without settup up parallelization
    for i = 1:length(cvf_params)
        % Print useful info about the model about to be fitted.
        fprintf('\n=================\nGroup %d of %d\n=================\n',...
                cvf_params(i).groupIdx, numGroups);
        if cvf_params(i).cvf == 0
            fprintf('\n===== Training on all data =====\n');
        else
            fprintf('\n===== Cross-validation fold %d of %d =====\n', ...
                cvf_params(i).cvf, numFolds);
        end
        fprintf('Number of training trials: %d\n', length(cvf_params(i).seqTrain));
        fprintf('Number of test trials: %d\n', length(cvf_params(i).seqTest));
        fprintf('Latent dimensionality: %d\n', cvf_params(i).xDim);
        fprintf('Observation dimensionality: %d\n', size(cvf_params(i).seqTrain(1).y,1));
        
        % Call the FA engine
        call_fa_engine(cvf_params(i).fname, cvf_params(i).seqTrain, ...
            cvf_params(i).seqTest, 'method', method, ...
            'xDim', cvf_params(i).xDim, 'cvf', cvf_params(i).cvf, ...
            'groupIdx', cvf_params(i).groupIdx, 'parallelize', parallelize, ...
            'binWidth', binWidth, 'rngSettings', rngSettings, extraOpts{:});        
    end
end