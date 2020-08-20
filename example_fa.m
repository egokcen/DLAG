% =========
% FA DEMO
% =========
%
% This demo shows how we can use FA to study multi-population data. It's 
% recommended to run this script section-by-section, rather than all at 
% once (or put a break point before Section 2, as it may take a long time).
%
% Section 1 demonstrates how FA can be used for exploratory data
% analysis.
%
%     Section 1a fits (to each group indepently) FA models with specified 
%     numbers of latent dimensions. Optional arguments are explicitly 
%     specified for the sake of demonstration.
%
%     Section 1b takes this model and performs basic inference of
%     latent trajectories. It also demonstrates how to orthonormalize and
%     order FA latents according to shared variance explained.
%
% Section 2 shows how to select the optimal dimensionality using
% cross-validation (for each group independently, as in Section 1).
%
% Author:
%     Evren Gokcen    egokcen@cmu.edu
%
% Last Revised:
%     16 Aug 2020

%% ==================
% 0a) Load demo data
% ======================

% Synthetic data generated from a pCCA model
dat_file = 'mat_sample/dlag_demo_data_synthetic';
fprintf('Reading from %s \n',dat_file);
load(dat_file);

%% =======================
% 0b) Set up parallelization
% ===========================

% If parallelize is true, all cross-validation folds will be analyzed in 
% parallel using Matlab's parfor construct. If you have access to multiple 
% cores, this provides significant speedup.
parallelize = false;
numWorkers = 2;      % Adjust this to your computer's specs

%% =====================
% 1a) Fitting a FA model
% ========================

% Let's explicitly define all of the optional arguments, for
% the sake of demonstration:
runIdx = 4;           % Results will be saved in baseDir/mat_results/runXXX/, where
                      % XXX is runIdx. Use a new runIdx for each dataset.
baseDir = '.';        % Base directory where results will be saved
overwriteExisting = true; % Control whether existing results files are overwritten
saveData = false;     % Set to true to save train and test data (not recommended)
binWidth = 20;        % Sample period / spike count bin width, in units of time (e.g., ms)
numFolds = 0;         % Number of cross-validation folds (0 means no cross-validation)
xDims = {4, 4};       % The number of latents for each group
yDims = [10 10];      % Number of observed features (neurons) in each group (area)
maxIters = 1e8;       % Limit the number of EM iterations (not recommended, in general)
randomSeed = 0;       % Seed the random number generator, for reproducibility
numGroups = length(yDims); % Number of observation groups in the data

fit_fa(runIdx, Ytrain, ...
       'baseDir', baseDir, ...
       'binWidth', binWidth, ...
       'numFolds', numFolds, ...
       'xDims', xDims, ...
       'yDims', yDims, ...
       'maxIters', maxIters, ...
       'parallelize', false, ...  % Only relevant if cross-validating
       'randomSeed', randomSeed, ...
       'overwriteExisting', overwriteExisting, ...
       'saveData', saveData);

%% =====================================
% 1b) Explore extracted latent trajectories
% ==========================================

% Convert data in seq format to format compatible with FA functions
Ytest_static = seq2pcca(Ytest, yDims, 'datafield', 'y'); % Borrowing code from pCCA pack

for groupIdx = 1:numGroups
    % Retrieve the fitted model of interest
    res = getModel_fa(runIdx, groupIdx, xDim, 'baseDir', baseDir);

    % Infer unordered latents
    [Xinferred, ~] = fa_estep(Ytest_static{groupIdx}, res.estParams);
    % Convert inferred latents back into seq format
    Ytest = segmentByTrial(Ytest, Xinferred.mean, 'xsm');
    % Plot unordered latents vs time
    plotEachDimVsTime(Ytest, 'xsm', res.binWidth, ...
                      'nPlotMax', 1, ...
                      'nCol', res.xDim, ...
                      'plotSingle', true, ...
                      'plotMean', true, ...
                      'units', 'ms');
                  
    % Orthonormalize and order latents according to shared variance
    % explained.
    xspec = sprintf('xorth%02d', res.xDim);
    [Xorth, Corth] = orthogonalize(Xinferred.mean, res.estParams.C);
    Ytest = segmentByTrial(Ytest, Xorth, xspec);
    plotEachDimVsTime(Ytest, xspec, res.binWidth, ...
                      'nPlotMax', 20, ...
                      'nCol', res.xDim, ...
                      'plotSingle', true, ...
                      'plotMean', true, ...
                      'units', 'ms');

    % Visualize the top three orthonormalized latents in 3D space
    plotTraj(Ytest, xspec, ...
             'dimsToPlot', 1:3, ...
             'nPlotMax', 1, ...
             'plotSingle', true, ...
             'plotMean', true);
end

%% ==========================
% 2) Cross-validate FA models
%  =============================

% Change other input arguments as appropriate
runIdx = 5;
numFolds = 4;
xDims = {0:yDims(1)-1, 0:yDims(2)-1}; % Sweep over these dimensionalities

fit_fa(runIdx, Ytrain, ...
       'baseDir', baseDir, ...
       'binWidth', binWidth, ...
       'numFolds', numFolds, ...
       'xDims', xDims, ...
       'yDims', yDims, ...
       'maxIters', maxIters, ...
       'parallelize', parallelize, ...
       'randomSeed', randomSeed, ...
       'numWorkers', numWorkers, ...
       'overwriteExisting', overwriteExisting, ...
       'saveData', saveData);

%% Get cross-validation results
[cvResults, bestModels] = getCrossValResults_fa(runIdx, 'baseDir', baseDir);

% Plot cross-validated performance vs estimated dimensionality
plotPerfvsDim_fa(cvResults, ...
                 'bestModels', bestModels);
               
% Find a conservative estimate of optimal dimensionality.
cutoffPC = 0.95;
d_shared = nan(1,numGroups);
for groupIdx = 1:numGroups
    d_shared(groupIdx) ...
        = findSharedDimCutoff_fa(cvResults{groupIdx}(bestModels(groupIdx)).estParams, ...
                                 cutoffPC, ...
                                 'plotSpec', true);
end
        