function result = fit_dlag(runIdx, dat, varargin)
%
% result = fit_dlag(runIdx, dat, ...)
% 
% Description: This function does most of the front-end work for fitting
%              the DLAG model. It prepares directories for saved results,
%              sets up cross-validation folds (along with parallelization),
%              and then calls the functions that actually estimate DLAG 
%              parameters.
%
%              NOTE: fit_dlag will not train models for which the number of
%                    within- and across-group dimensions adds up to more
%                    than the number of features (neurons) within a group
%                    (area). You can input such models in xDims_across and
%                    xDims_within, but they will be ignored here.
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
%     method      -- string; method to be used (currently one supported):
%                    'dlag'. (default: 'dlag')
%     binWidth    -- float; bin width or sample period, in units of time
%                    (e.g., ms) (default: 20)
%     numFolds    -- int; number of cross-validation folds (default: 0).
%                    0 indicates no cross-validation, i.e. train on all trials.
%     yDims       -- (1 x numGroups) array; Specify the number of features 
%                    (neurons) in each group (area). Elements in yDims
%                    should match the format of dat.
%     xDims_across -- (1 x numDims) array; across-group state 
%                     dimensionalities to be modeled. For example, to fit
%                     3 separate models, each with 1, 2, and 3 across-group
%                     dimensions, input 1:3 or [1 2 3]. (default: [3])
%     xDims_within -- (1 x numGroups) cell array; each element is a vector of
%                     within-group state dimensionalities to be modeled, of
%                     the same format as xDims_across. For example, {1:3,
%                     1:3} specifies desired within-group state
%                     dimensionalities for data with two groups.
%                     If xDims_within contains only one element, but there
%                     are multiple groups (as indicated by yDims,
%                     xDims_across), then all groups will be assigned the
%                     same within-group dimensionalities specified in that
%                     vector. (default: {[1]})
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
%               method      -- string; method that was used (for now, just
%                              'dlag')
%               rngSettings -- structure with the random number generator 
%                              settings used during run time. Includes 
%                              fields 'Type', 'Seed', and 'State'.
%               init_method -- string; method used to initialize DLAG
%                              ('pCCA', 'params')
%               xDim_across -- int; number of across-group dimensions
%               xDim_within -- (1 x numGroups) array; number of within-
%                              group dimensions for each group
%               yDims       -- (1 x numGroups) array; observation
%                              dimensionalities of each group
%               rGroups     -- (1 x 2) array; the indexes of two groups 
%                              used to measure generalization performance
%                              via regression
%               binWidth    -- float; bin width or sample period, in units 
%                              of time (e.g., ms)
%               startTau    -- float; value used to initialize GP 
%                              timescales, in same units of time as 
%                              binWidth (not relevant if init_method is
%                              'params')
%               startEps    -- float; value used to initialize GP noise
%                              variances (not relevant if init_method is 
%                              'params')
%               startDelay  -- (1 x numGroups) array; Values used to 
%                              initialize delays between across-area  
%                              latents and groups. (NOT the same format as  
%                              the DelayMatrix DLAG parameter.) Entries in 
%                              units of time. (not relevant if init_method 
%                              is 'params')
%               covType     -- string; Type of GP covariance kernel used 
%                              (e.g., 'rbf')
%               parallelize -- logical; indicates whether or not 
%                              parallelization was used.
%               saveData    -- logical; indicates whether or not train and 
%                              test data is saved in this file.
%               startParams -- structure containing parameter values
%                              at which DLAG EM procedure was initialized.
%               estParams   -- structure containing final DLAG parameter
%                              estimates, after model fitting was completed.
%               iterTime    -- (1 x numIters) array; iterTime(i) contains
%                              the amount of clock time, in sec, EM 
%                              iteration i took to complete.
%               D           -- (1 x numIters) cell array; the estimated 
%                              DLAG delay matrix after each EM iteration.
%               gams_across -- (1 x numIters) cell arry; estimated 
%                              gamma_across after each EM iteration.
%               gams_within -- (1 x numGroups) cell arry;
%                              gams_within(i) -- (1 x numIters) cell array;
%                              estimated gamma_within for group i after 
%                              each EM iteration.
%               LLcut       -- (1 x numIters) array; data log likelihood  
%                              after each EM iteration (where training data 
%                              was potentially 'cut' into trials of equal 
%                              length)
%                              Note: Entries will be NaN, on iterations
%                                    was LL was not computed, to save time.
%               LLtrain     -- float; final data log likelihood evaluated 
%                              on uncut train data
%               cvf         -- int; indicates which cross-validation fold 
%                              this model was fit to. (0 if fit to all 
%                              data)
%               seqTrain    -- training data structure of the same format 
%                              as dat. Also includes latent trajectories 
%                              inferred using estParams.
%                              NOTE: Only included if saveData is true. 
%               seqTrainCut -- data structure of the same format as
%                              seqTrain, but for data where trials were
%                              cut to be the same length.
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
%     MSEorth -- (length(mList) x 2) array; mean-squared error in each 
%                direction for reduced DLAG predictions, for pairwise
%                regression performed on test data
%     R2      -- (1 x 2) array; R^2 in each pairwise direction, for
%                pairwise regression performed on test data
%     R2orth  -- (length(mList) x 2) array; R^2 error in each 
%                direction for reduced DLAG predictions, for pairwise
%                regression  performed on test data
%     seqTrainCut -- test data structure of the same format as seqTrain.
%                    NOTE: Only included if saveData is true.
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     11 Mar 2020 -- Initial full revision.
%     07 May 2020 -- Added option to seed random number generator.
%     13 May 2020 -- Error will be thrown if yDims specified incorrectly.
%                    Removed hasSpikesBool functionality, which removed
%                    inactive units based on the training set.
%     14 May 2020 -- Removed datFormat option. Spiking data can be
%                    preprocessed separately, if desired, with getSeq.m

% Specify defaults for optional arguments
baseDir              = '.';
method               = 'dlag';
binWidth             = 20;
numFolds             = 0;
yDims                = [];
xDims_across         = [3];
xDims_within         = {[1]};
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

% Construct hyperparameter grid. Note that models will not be allowed for
% which latent dimensions outnumber observations. (These models will have 
% issues with numerical conditioning.)
% Dimensions len(xDims_across) x len(xDims_within{1}) x len(xDims_within{2}) ...
xDims_within_grid = cell(1,length(xDims_within));
[xDims_across_grid, xDims_within_grid{:}] = ndgrid(xDims_across,xDims_within{:});
xDims_grid = [xDims_across_grid(:)];
for groupIdx = 1:numGroups
    if length(xDims_within) == 1
        % If user input only one set of private dimensions, then we'll fix
        % the private dimensionalities of all groups to be the same.
        xDims_grid = [xDims_grid xDims_within_grid{1}(:)];
    else
        % Otherwise, do a full grid search
        xDims_grid = [xDims_grid xDims_within_grid{groupIdx}(:)];
    end
end
for paramIdx = 1:size(xDims_grid,1)
    xDims(paramIdx).xDim_across = xDims_grid(paramIdx,1);
    xDims(paramIdx).xDim_within = xDims_grid(paramIdx,2:end);
end

% Set cross-validation folds
for cvf = cvf_list 
    val = false(1, N);
    if cvf > 0
        val = (val_indices == cvf);
    end
    train = ~val;

    seqTrain      = dat(train);
    seqTest       = dat(val);

    % Check if training data covariance is full rank
    yAll = [seqTrain.y];
    yDim  = size(yAll, 1);
    if rank(cov(yAll')) < yDim
        fprintf('ERROR: Observation covariance matrix is rank deficient.\n');
        fprintf('Possible causes: repeated units, not enough observations.\n');
        fprintf('Exiting...\n');
        return
    end

    cv_data(cvf+1).seqTrain = seqTrain;
    cv_data(cvf+1).seqTest = seqTest;
end

% We assume within- and across-group dimensions are linearly 
% independent, so don't train models where a group has an 
% overcomplete basis (i.e., more latents than observations).
i = 1;    
for paramIdx = 1:length(xDims)
    % Flag whether or not we should skip this set of hyperparameters
    overcomplete = false; 
    for groupIdx = 1:numGroups
        if xDims(paramIdx).xDim_across + xDims(paramIdx).xDim_within(groupIdx) > yDims(groupIdx)
            overcomplete = true;
        end
    end
    
    % If overcomplete, then we skip this set of hyperparameters
    if ~overcomplete
    
        for cvf = cvf_list
            % Specify filename where results will be saved
            fname = sprintf('%s/%s_nGroups%02d_xDimA%02d_xDimW', runDir, method, numGroups, xDims(paramIdx).xDim_across);
            for groupIdx = 1:numGroups
                fname = sprintf('%s_%02d', fname, xDims(paramIdx).xDim_within(groupIdx)); 
            end
            if cvf > 0
                fname = sprintf('%s_cv%02d', fname, cvf);
            end
            cvf_params(i).fname = fname;
            cvf_params(i).xDim_across = xDims(paramIdx).xDim_across;
            cvf_params(i).xDim_within = xDims(paramIdx).xDim_within;
            cvf_params(i).cvf = cvf;
            i = i+1;
        end
    
    end
end

% Set up parallelization, if desired
if parallelize
    StartParPool(numWorkers);  % Helper function to set up parfor construct
    parfor i = 1:length(cvf_params)
        call_dlag_engine(cvf_params(i).fname, cv_data(cvf_params(i).cvf+1).seqTrain, ...
            cv_data(cvf_params(i).cvf+1).seqTest, 'method', method, ...
            'xDim_across', cvf_params(i).xDim_across, 'xDim_within', cvf_params(i).xDim_within,...
            'yDims', yDims, 'cvf', cvf_params(i).cvf, 'parallelize', parallelize, ...
            'binWidth', binWidth, 'rngSettings', rngSettings, extraOpts{:});                
    end
else % Otherwise, continue without settup up parallelization
    for i = 1:length(cvf_params)
        if cvf_params(i).cvf == 0
            fprintf('\n===== Training on all data =====\n');
        else
            fprintf('\n===== Cross-validation fold %d of %d =====\n', ...
                cvf_params(i).cvf, numFolds);
        end
        call_dlag_engine(cvf_params(i).fname, cv_data(cvf_params(i).cvf+1).seqTrain, ...
            cv_data(cvf_params(i).cvf+1).seqTest, 'method', method, ...
            'xDim_across', cvf_params(i).xDim_across, 'xDim_within', cvf_params(i).xDim_within, ...
            'yDims', yDims, 'cvf', cvf_params(i).cvf, 'parallelize', parallelize, ...
            'binWidth', binWidth, 'rngSettings', rngSettings, extraOpts{:});        
    end
end

% Results are saved to a file.
% Here, return those results (from the file) only if not doing
% cross-validation
result = {};
if (nargout == 1) && (numFolds == 0) && exist([cvf_params(1).fname '.mat'], 'file')
    result = load(cvf_params(1).fname);
end