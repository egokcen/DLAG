% =========
% pCCA DEMO 
% =========
%
% This demo shows how we can extract latent variables from multi-population
% data with static methods like pCCA. It's recommended to run this script
% section-by-section, rather than all at once (or put a break point before
% Section 2, as it may take a long time).
%
% Section 1 demsonstrates how pCCA can be used for exploratory data 
% analysis.
%
%     Section 1a fits a pCCA model with a specified number of across-group
%     latent dimensions. Optional arguments are explicitly specified for 
%     the sake of demonstration.
%
%     Section 1b takes this model and performs basic inference of 
%     latent trajectories. It also demonstrates how to orthonormalize and
%     order pCCA latents according to across-group predictive power.
%
% Section 2 shows how to select the optimal dimensionality using 
% cross-validation. 
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Last Revised: 
%     11 Apr 2020

%% ==================
% 0) Load demo data
% ======================

% Synthetic data generated from a pCCA model
dat_file = 'mat_sample/pcca_demo_data_synthetic';
fprintf('Reading from %s \n',dat_file);
load(dat_file);

%% =====================
% 1a) Fitting a pCCA model
% ========================

% Let's explicitly define all of the optional arguments, for 
% the sake of demonstration:
runIdx = 4;           % Results will be saved in baseDir/mat_results/runXXX/, where 
                      % XXX is runIdx. Use a new runIdx for each dataset.
datFormat = 'seq';    % Analyzing continuous valued data (e.g., spike counts)
baseDir = '.';        % Base directory where results will be saved
overwriteExisting = true; % Control whether existing results files are overwritten
saveData = false;     % Set to true to save train and test data (not recommended)
method = 'pcca';      % Specify static method used (only pCCA, for now)
binWidth = 20;        % Sample period / spike count bin width, in units of time (e.g., ms)
numFolds = 0;         % Number of cross-validation folds (0 means no cross-validation)
xDim = 4;             % This number of across-group latents matches the synthetic ground truth
yDims = [10 10];      % Number of observed features (neurons) in each group (area)
rGroups = [1 2];      % For performance evaluation, we can regress Area 2's activity with Area 1
maxIters = 1e8;       % Limit the number of EM iterations (not recommended, in general)
parallelize = false;  % Only relevant if cross-validating

fit_pcca(runIdx, Ytrain, ...
         'datFormat', datFormat, ...
         'baseDir', baseDir, ...
         'method', method, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'xDim', xDim, ...
         'yDims', yDims, ...
         'rGroups', rGroups,...
         'maxIters', maxIters, ...
         'parallelize', parallelize, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);

%% =====================================
% 1b) Explore extracted latent trajectories
% ==========================================

% Retrieve the fitted model of interest
res = getModel_pcca(runIdx, xDim, 'baseDir', baseDir);

% Convert data in seq format to format compatible with pCCA functions
Ys_test = seq2pcca(Ytest.seq, res.yDims, 'datafield', 'y');
% Infer unordered latents
[Xinferred, ~] = pcca_estep(Ys_test, res.estParams);
% Convert inferred latents back into seq format
Ytest.seq = segmentByTrial(Ytest.seq, Xinferred.mean, 'xsm');
% Plot unordered latents vs time
plotEachDimVsTime(Ytest.seq, 'xsm', res.binWidth, ... 
                  'nPlotMax', 1, ...
                  'nCol', res.xDim, ...
                  'plotSingle', true, ...
                  'plotMean', true, ...
                  'units', 'ms');
              
% Orthonormalize and order across-area latents according to predictive
% power.
xspec = sprintf('xorth%02d', res.xDim);
[Xorth, Corth] = predictiveProjection_pcca(Ys_test, ...
                                           res.estParams, ... 
                                           res.xDim, ...
                                           res.rGroups);
Ytest.seq = segmentByTrial(Ytest.seq, Xorth, xspec);
plotEachDimVsTime(Ytest.seq, xspec, res.binWidth, ...
                  'nPlotMax', 1, ...
                  'nCol', res.xDim, ...
                  'plotSingle', true, ...
                  'plotMean', true, ...
                  'units', 'ms');

% Visualize the top three orthonormalized latents in 3D space
plotTraj(Ytest.seq, xspec, ...
         'dimsToPlot', 1:3, ...
         'nPlotMax', 1, ...
         'plotSingle', true, ...
         'plotMean', true);

%% ==========================
% 2) Cross-validate pCCA models
%  =============================

% If parallelize is true, all folds will be run in parallel using Matlab's
% parfor construct. If you have access to multiple cores, this provides
% significant speedup.
parallelize = false;
numWorkers = 2;      % Adjust this to your computer's specs

% Change other input arguments as appropriate
runIdx = 5;
numFolds = 4;
xDim = 0:min(yDims); % Sweep over these dimensionalities

fit_pcca(runIdx, Ytrain, ...
         'datFormat', datFormat, ...
         'baseDir', baseDir, ...
         'method', method, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'xDim', xDim, ...
         'yDims', yDims, ...
         'rGroups', rGroups,...
         'maxIters', maxIters, ...
         'parallelize', parallelize, ...
         'numWorkers', numWorkers, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);
      
%% Get cross-validation results
[cvResults, bestModel] = getCrossValResults_pcca(runIdx, 'baseDir', baseDir);

% Plot cross-validated performance vs estimated dimensionality
plotPerfvsDim_pcca(cvResults, ...
                   'bestModel', bestModel, ...
                   'plotLL', true, ...
                   'plotR2', false, ...
                   'plotMSE', false);
