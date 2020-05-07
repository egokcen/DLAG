function result = fit_pcca(runIdx, dat, varargin)
%
% result = fit_pcca(runIdx, dat,...)
%
% Description: This function does most of the front-end work for fitting
%              the probabilistic canonical correlation analysis(pCCA).
%              It prepares data, sets up cross-validation folds (along with
%              parallelization), and then calls the functions that actually
%              estimate model parameters.
%              Data can be input in one of two formats: as a 0/1 matrix of 
%              spiking activity or a sequence of continuous values (see 
%              datFormat).
%
% Arguments: 
%
%     Required:
%
%     runIdx    --  int; results files will be saved in 
%                   baseDir/mat_results/runXXX, where XXX is runIdx.
%                   baseDir can be specified by the user (see below)
%     dat       --  (1 x N) structure whose nth entry (corresponding to the
%                   nth experimental trial) has fields that depend on the 
%                   value of datFormat. In all cases, features (neurons)
%                   are assumed to be ordered according to which group
%                   (area) they belong to. For example, if data has
%                   features belonging to two groups, then each observation
%                   should be organized as follows: 
%                        [      |
%                         Group 1 Features
%                               |
%                               |
%                         Group 2 Features
%                               |         ]
%                   if datFormat is 'spikes', nth entry has fields
%                       trialId -- unique trial identifier                       
%                       spikes  -- (yDim x T) array; 0/1 matrix of the 
%                                  raw spiking activity across all 
%                                  neurons.  Each row corresponds to a 
%                                  neuron. Each column corresponds to a
%                                  1 msec timestep.
%                   if datFormat is 'seq', nth entry has fields
%                       trialId      -- unique trial identifier  
%                       T (1 x 1)    -- number of timesteps
%                       y (yDim x T) -- continuous valued data 
%                                       (Eg: binned spike counts)
%     yDims     -- (1 x numGroups) array; Specify the number features 
%                  (neurons) in each group (area). Elements in yDims
%                  should match the format of dat. 
%                  NOTE: This is a named argument, for convenience of
%                        implementation, but the code will break if it is
%                        not specified.
%
%     Optional: (Many of these are technically defined as optional, but
%               should really be specified).
%
%     datFormat   -- string; specifies format of input data:
%                    'spikes' or 'seq'. (default: 'seq')
%     baseDir     -- string; specifies directory in which to store
%                    mat_results. (default: '.', i.e., current directory)
%     method      -- string; method to be used (right now, just 'pcca') 
%                    (default: 'pcca')
%     binWidth    -- float; bin width or sample period, in units of time
%                    (e.g., ms) (default: 20)
%     numFolds    -- int; number of cross-validation folds (default: 0).
%                    0 indicates no cross-validation, i.e. train on all trials.
%     xDims       -- (1 x numDims) array; across-group state 
%                     dimensionalities to be modeled. For example, to fit
%                     3 separate models, each with 1, 2, and 3 across-group
%                     dimensions, input 1:3 or [1 2 3]. (default: [3])
%     rGroups     -- (1 x 2) array; the indexes of two groups 
%                    used to measure generalization performance
%                    via regression (default: [1 2])
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
%     result -- structure containing all variables saved in 
%               ./mat_results/runXXX/ if 'numFolds' is 0.  
%               Else, the structure is empty. These variables should help
%               with reproducibility, and include:
%               
%               fname       -- string; path to file where variables are 
%                              stored
%               method      -- string; method that was used ('pcca')
%               rngSettings -- structure with the random number generator 
%                              settings used during run time. Includes 
%                              fields 'Type', 'Seed', and 'State'.
%               xDim        -- int; number of across-group dimensions
%               yDims       -- (1 x numGroups) array; observation
%                              dimensionalities of each group
%               rGroups     -- (1 x 2) array; the indexes of two groups 
%                              used to measure generalization performance
%                              via regression
%               binWidth    -- float; bin width or sample period, in units 
%                              of time (e.g., ms)
%               parallelize -- logical; indicates whether or not 
%                              parallelization was used.
%               saveData    -- logical; indicates whether or not train and 
%                              test data is saved in this file.
%               startParams -- structure containing parameter values
%                              at which pCCA EM procedure was initialized.
%               estParams   -- structure containing final static model 
%                              parameter estimates, after model fitting was
%                              completed.
%               LL          -- (1 x numIters) array; data log likelihood  
%                              after each pCCA EM iteration
%               LLtrain     -- float; final data log likelihood evaluated 
%                              on train data
%               cvf         -- int; indicates which cross-validation fold 
%                              this model was fit to. (0 if fit to all 
%                              data)
%               hasSpikesBool -- (yDim x 1) logical array; indicates which
%                                neurons had non-zero activity.
%               seqTrain    -- training data structure of the same format 
%                              as dat (with data formatted as in the 'seq' 
%                              format option). Also includes latent 
%                              trajectories inferred using estParams.
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
%     MSE     -- (1 x 2) array; mean-squared error in each pairwise 
%                direction, for pairwise regression performed on test data
%     R2      -- (1 x 2) array; R^2 in each pairwise direction, for
%                pairwise regression performed on test data
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     11 Mar 2020 -- Initial full revision.
%     07 May 2020 -- Added option to seed random number generator.

% Specify defaults for optional arguments
datFormat            = 'seq';
baseDir              = '.';
method               = 'pcca';
binWidth             = 20;
numFolds             = 0;
yDims                = [];
xDims                = [3];
rGroups              = [1 2];
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

% If input is in format 1 (raw spikes), obtain binned spike counts
switch(datFormat)
    case 'spikes'
        seq  = getSeq(dat, binWidth, extraOpts{:});
    case 'seq'
        seq  = dat.seq;
end

if isempty(seq)
    fprintf('Error: No valid trials.  Exiting.\n');
    return;
end

N    = length(seq);      % Number of trials
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
for xDim = xDims
    for cvf = cvf_list
        % Set cross-validation folds
        val = false(1, N);
        if cvf > 0
            val = (val_indices == cvf);
        end
        train = ~val;
        
        seqTrain      = seq(train);
        seqTest       = seq(val);
        
        % Remove inactive units based on training set
        hasSpikesBool = (mean([seqTrain.y], 2) ~= 0); %1:size(seqTrain(1).y,1); %
        
        for n = 1:length(seqTrain)
            seqTrain(n).y = seqTrain(n).y(hasSpikesBool,:);
        end
        
        for n = 1:length(seqTest)
            seqTest(n).y = seqTest(n).y(hasSpikesBool,:);
        end
        
        % Check if training data covariance is full rank
        yAll = [seqTrain.y];
        yDim  = size(yAll, 1);
        if rank(cov(yAll')) < yDim
            fprintf('ERROR: Observation covariance matrix is rank deficient.\n');
            fprintf('Possible causes: repeated units, not enough observations.\n');
            fprintf('Exiting...\n');
            return
        end
        
        cvf_params(i).seqTrain = seqTrain;
        cvf_params(i).seqTest = seqTest;
        cvf_params(i).hasSpikesBool = hasSpikesBool;
        
        % Specify filename where results will be saved
        fname = sprintf('%s/%s_xDim%02d', runDir, method, xDim);
        if cvf > 0
            fname = sprintf('%s_cv%02d', fname, cvf);
        end
        cvf_params(i).fname = fname;
        cvf_params(i).xDim = xDim;
        cvf_params(i).cvf = cvf;
        i = i+1;
    end
end

% Set up parallelization, if desired
if parallelize
    StartParPool(numWorkers);  % Helper function to set up parfor construct
    parfor i = 1:length(cvf_params)
        call_pcca_engine(cvf_params(i).fname, cvf_params(i).seqTrain, ...
            cvf_params(i).seqTest, 'method', method, 'yDims', yDims, ...
            'xDim', cvf_params(i).xDim, 'cvf', cvf_params(i).cvf, ...
            'rGroups', rGroups, 'parallelize', parallelize, ...
            'hasSpikesBool', cvf_params(i).hasSpikesBool, 'binWidth', binWidth,...
            'rngSettings', rngSettings, extraOpts{:});                
    end
else % Otherwise, continue without settup up parallelization
    for i = 1:length(cvf_params)
        if cvf_params(i).cvf == 0
            fprintf('\n===== Training on all data =====\n');
        else
            fprintf('\n===== Cross-validation fold %d of %d =====\n', ...
                cvf_params(i).cvf, numFolds);
        end
        call_pcca_engine(cvf_params(i).fname, cvf_params(i).seqTrain, ...
            cvf_params(i).seqTest, 'method', method, 'yDims', yDims, ...
            'xDim', cvf_params(i).xDim, 'cvf', cvf_params(i).cvf, ...
            'rGroups', rGroups, 'parallelize', parallelize, ...
            'hasSpikesBool', cvf_params(i).hasSpikesBool, 'binWidth', binWidth,...
            'rngSettings', rngSettings, extraOpts{:});        
    end
end

% Results are saved to a file.
% Here, return those results (from the file) only if not doing
% cross-validation
result = {};
if (nargout == 1) && (numFolds == 0) && exist([cvf_params(1).fname '.mat'], 'file')
    result = load(cvf_params(1).fname);
end