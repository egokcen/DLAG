% =========
% DLAG DEMO 
% ========= 
%
% This demo shows how we can extract latent variables from multi-population
% data with DLAG (Gokcen et al., 2020). It's recommended to run this script
% section-by-section, rather than all at once (or put a break point before
% Sections 2 and 3, as they may take a long time, depending on your use
% of parallelization).
%
% Section 1 demonstrates how DLAG can be used for exploratory data 
% analysis.
%
%     Section 1a fits a DLAG model with a specified number of within-
%     and across-group latent dimensions. Optional arguments are
%     explicitly specified for the sake of demonstration.
%
%     Section 1b takes this model and explores the latent GP timescales and
%     delays. It performs basic inference of within- and across-group 
%     latent trajectories. One can compare estimated parameters and
%     trajectories to the ground truth that underlies the demo synthetic
%     data.
%
%     Section 1c demonstrates how to orthonormalize and order latents 
%     according to various objectives, such as predictive power across
%     groups or shared variance explained within a group.
%
%     Section 1d demonstrates how to denoise observations using a DLAG
%     model. One can compare raw observations to their denoised
%     counterparts.
%
% Section 2 shows how to select optimal DLAG dimensionalities using 
% a streamlined cross-validation approach. 
%
%     Section 2a uses cross-validation to first determine the optimal 
%     across-group dimensionality according to a pCCA model (see 
%     example_pcca.m for a detailed demo of pCCA). Applying pCCA or any 
%     related static method first is recommended. It is much faster than 
%     DLAG, and if no significant interaction is found (i.e., pCCA 
%     cross-validation returns 0 across-group dimensions), then proceed no
%     further. The pCCA across-group dimensionality estimate can be used in
%     tandem with the FA dimensionality estimates found in Section 2b for 
%     initial exploratory analysis with DLAG, prior to more careful
%     cross-validation--demonstrated in Section 2c.
%
%     Section 2b establishes an upper bound on the total dimensionality of 
%     each group (within + across) by applying FA to each group 
%     independently. See example_fa.m for a detailed demo of FA on
%     multi-population data.
%
%     Section 2c determines the optimal DLAG across- and within-group 
%     dimensionalities using a streamlined cross-validation approach. 
%     The search space is constrained to models such that the across-
%     and within-group dimensionalities for each group add up to the upper
%     bounds established by FA in Section 2b. Leave-population-out
%     predictive performance is used to choose the model with the fewest
%     across-group dimensions that still exhibits close-to-optimal 
%     performance. 
%
% Section 3 demonstrates several post-selection inference procedures.
% After selecting the optimal model via cross-validation, these procedures
% can elucidate the prominence of individual latent variables and
% the uncertainty in parameter estimates.
%
%     Section 3a analyzes the spectrum of the DLAG loading matrices, to 
%     determine the number of dimensions needed to capture a certain 
%     proportion of the shared variance, among all groups jointly, 
%     across groups, and within each group individually.
%
%     Section 3b evaluates the bootstrapped prominence of individual latent
%     variables--across groups jointly and within each group individually.
%
%     Section 3c evaluates how significantly each across-group set of 
%     delays deviates from 0, using bootstrapped samples.
%
%     Section 3d constructs bootstrapped confidence intervals for latent
%     delays and timescales. This section involves re-fitting DLAG models
%     to bootstrapped samples, so its runtime is similar to that of 
%     Section 2, depending on how much parallelization is used.
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
%
% Last Revised: 
%     18 Aug 2020

%% ================
% 0a) Load demo data 
% ===================

% Synthetic data generated from a DLAG model
dat_file = 'mat_sample/dlag_demo_data_synthetic';
fprintf('Reading from %s \n',dat_file);
load(dat_file);

%% =======================
% 0b) Set up parallelization
% ===========================

% If parallelize is true, all cross-validation folds and bootstrap samples
% will be analyzed in parallel using Matlab's parfor construct. 
% If you have access to multiple cores, this provides significant speedup.
parallelize = false;
numWorkers = 2;      % Adjust this to your computer's specs

%% =====================
% 1a) Fitting a DLAG model
% ========================

% Let's explicitly define all of the optional arguments, for 
% the sake of demonstration:
runIdx = 1;           % Results will be saved in baseDir/mat_results/runXXX/,  
                      % where XXX is runIdx. Use a new runIdx for each dataset.
baseDir = '.';        % Base directory where results will be saved
overwriteExisting = false; % Control whether existing results files are overwritten
saveData = false;     % Set to true to save train and test data (not recommended)
method = 'dlag';      % For now this is the only option, but that may change in the near future
binWidth = 20;        % Sample period / spike count bin width, in units of time (e.g., ms)
numFolds = 0;         % Number of cross-validation folds (0 means no cross-validation)
xDims_across = 4;     % This number of across-group latents matches the synthetic ground truth
xDims_within = {2, 2}; % These numbers match the within-group latents in the synthetic ground truth
yDims = [10 10];      % Number of observed features (neurons) in each group (area)
rGroups = [1 2];      % For performance evaluation, we can regress group 2's activity with group 1
startTau = 2*binWidth;% Initial timescale, in the same units of time as binWidth
segLength = 25;       % Largest trial segment length, in no. of time points
init_method = 'static'; % Initialize DLAG with fitted pCCA parameters
learnDelays = true;   % Set to false if you want to fix delays at their initial value
maxIters = 5e3;       % Limit the number of EM iterations (not recommended, in general)
freqLL = 10;          % Check for data log-likelihood convergence every freqLL EM iterations
freqParam = 100;      % Store intermediate delay and timescale estimates every freqParam EM iterations
minVarFrac = 0.01;    % Private noise variances will not be allowed to go below this value
verbose = true;       % Toggle printed progress updates
randomSeed = 0;       % Seed the random number generator, for reproducibility

fit_dlag(runIdx, Ytrain, ...
         'baseDir', baseDir, ...
         'method', method, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'xDims_across', xDims_across, ...
         'xDims_within', xDims_within, ...
         'yDims', yDims, ...
         'rGroups', rGroups,...
         'startTau', startTau, ...
         'segLength', segLength, ...
         'init_method', init_method, ...
         'learnDelays', learnDelays, ...
         'maxIters', maxIters, ...
         'freqLL', freqLL, ...
         'freqParam', freqParam, ...
         'minVarFrac', minVarFrac, ...
         'parallelize', false, ... % Only relevant for cross-validation
         'verbose', verbose, ...
         'randomSeed', randomSeed, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);

%% =========================================================
% 1b) Explore extracted GP parameters and compare to ground truth
% ================================================================

% Retrieve the fitted model of interest
xDim_across = 4;
xDim_within = [2 2];
res = getModel_dlag(runIdx, xDim_across, xDim_within, ...
                    'baseDir', baseDir);

% Plot training progress of various quantities. These plots can help with
% troubleshooting, if necessary.
plotFittingProgress(res, ...
                    'freqLL', freqLL, ...
                    'freqParam', freqParam, ...
                    'units', 'ms');

% Plot estimated within-group GP timescales
plotGPparams_dlag(res.estParams, res.binWidth, res.rGroups, ...
                  'plotAcross', false, ...
                  'plotWithin', true, ...
                  'units', 'ms');
              
% Plot ground truth within-group GP timescales. 
% Note that latent variables are not, in general, ordered. So don't try to
% match estimated latent 1 to ground truth latent 1, and so on.
plotGPparams_dlag(trueParams, res.binWidth, res.rGroups, ...
                  'plotAcross', false, ...
                  'plotWithin', true, ...
                  'units', 'ms');

% Plot estimated and ground truth delays and across-group GP timescales
% together on the same plot. For these scatterplots, it's more
% straightforward to match ground truth latents to corresponding
% estimates.
plotGPparams_withGT_dlag(res.estParams, trueParams, res.binWidth,...
                         res.rGroups, 'units', 'ms');

% Plot unordered latents but maintain delay labels
[seqTrain, ~] = exactInferenceWithLL_dlag(Ytrain, res.estParams);
plotEachDimVsTime(seqTrain, 'xsm', res.binWidth, ...
                  'nPlotMax', 1, ...
                  'plotSingle', true, ...
                  'plotMean', true, ...
                  'units', 'ms', ...
                  'dlagFormat', true, ...
                  'params', res.estParams);

% Plot ground truth latents, in the same format as above.
% Note that latent variables are not, in general, ordered. So don't try to
% match estimated latent 1 to ground truth latent 1, and so on.
plotEachDimVsTime(Xtrain, 'xgt', res.binWidth, ...
                  'nPlotMax', 1, ...
                  'plotSingle', true, ...
                  'plotMean', true, ...
                  'units', 'ms', ...
                  'dlagFormat', true, ...
                  'params', trueParams);

%% ===============================================================
% 1c) Orthonormalize latent trajectories according to various objectives
% =======================================================================   

% Orthonormalize and order across-group latents according to predictive
% power. This loses information about individual delay/timescale labels.
[seqTrain, Cpred] = predictiveProjection_dlag(seqTrain, ...
                                              res.estParams, ...
                                              res.xDim_across, ...
                                              res.rGroups);

% Name of orthonormalized latents in seqTest  
xspec = sprintf('xorth_pred%02d', res.xDim_across);                                      
plotEachDimVsTime(seqTrain, xspec, res.binWidth, ...
                  'nPlotMax', 20, ...
                  'nCol', res.xDim_across, ...
                  'plotSingle', true, ...
                  'plotMean', true, ...
                  'units', 'ms');
              
% Visualize the top three orthonormalized latents in 3D space
plotTraj(seqTrain, xspec, ...
         'dimsToPlot', 1:3, ...
         'nPlotMax', 1, ...
         'plotSingle', true, ...
         'plotMean', true);
     
% Alternatively, orthonormalize and order latents according to overall 
% shared variance within each group. One can also focus on shared variance
% due only to across-group latents or within-group latents.
includeAcross = true;
includeWithin = true;
[seqTrain, Corth] = orthonormalizeWithinGroups(seqTrain, res.estParams, ...
                                              'includeAcross', includeAcross, ...
                                              'includeWithin', includeWithin);
xspec = 'xorth'; % This name changes depending on includeAcross/includeWithin
if includeAcross
    xspec = sprintf('%s_across', xspec);
end
if includeWithin
    xspec = sprintf('%s_within', xspec);
end                                 
plotEachDimVsTime(seqTrain, xspec, res.binWidth, ...
                  'nPlotMax', 20, ...
                  'nCol', res.xDim_across + res.xDim_within(1), ...
                  'plotSingle', true, ...
                  'plotMean', true, ...
                  'units', 'ms');                       

%% =======================================
% 1d) Denoise observations using a DLAG model
% ============================================

% Denoise observations
[seqTrain, ~, ~, ~, ~, ~] = denoise_dlag(seqTrain, res.estParams);

% Compare PSTHs of raw observations to PSTHs of denoised observations
psth_raw = get_psth(seqTrain, 'spec', 'y');
spec = sprintf('yDenoisedOrth%02d', sum(xDim_across + xDim_within));
psth_denoised = get_psth(seqTrain, 'spec', spec);
plotSeqRaster(psth_raw, res.binWidth, 'units', 'ms');
plotSeqRaster(psth_denoised, res.binWidth, 'units', 'ms');

%% =========================================================
% 2a) Cross-validate pCCA models to make sure there's significant
%     across-group interaction, and get an initial estimate of 
%     across-group dimensionality for exploratory analysis with
%     DLAG, if desired.
%  ===============================================================

% Change other input arguments as appropriate
runIdx = 2;
numFolds = 4;
xDim = 0:min(yDims)-1; % Sweep over these dimensionalities

fit_pcca(runIdx, Ytrain, ...
         'baseDir', baseDir, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'xDim', xDim, ...
         'yDims', yDims, ...
         'rGroups', rGroups,...
         'parallelize', parallelize, ...
         'numWorkers', numWorkers, ...
         'randomSeed', randomSeed, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);
      
%% Inspect full cross-validation results
[cvResults, bestModel] = getCrossValResults_pcca(runIdx, 'baseDir', baseDir);

% Plot cross-validated performance vs estimated dimensionality.
% Data log-likelihood is the standard performance metric. pCCA can also be
% used for regression between groups. To see cross-validated performance of
% pCCA as a regression method, set plotR2 and/or plotMSE to true.
plotPerfvsDim_pcca(cvResults, ...
                   'bestModel', bestModel, ...
                   'plotLL', true, ...
                   'plotR2', false, ...
                   'plotMSE', false);
               
% Find a conservative estimate of optimal across-group dimensionality.
cutoffPC = 0.95;
xDim_across_pcca = findSharedDimCutoff_pcca(cvResults(bestModel).estParams, ...
                                            cutoffPC, ...
                                            'plotSpec', true);
                                        
%% ===========================================================
% 2b) Cross-validate FA models to establish an upper bound on total
%     dimensionality (across+within) in each group.
%  =================================================================

% Change other input arguments as appropriate
runIdx = 2;
numFolds = 4;
xDims = {0:yDims(1)-1, 0:yDims(2)-1}; % Sweep over these dimensionalities

fit_fa(runIdx, Ytrain, ...
       'baseDir', baseDir, ...
       'binWidth', binWidth, ...
       'numFolds', numFolds, ...
       'xDims', xDims, ...
       'yDims', yDims, ...
       'parallelize', parallelize, ...
       'randomSeed', randomSeed, ...
       'numWorkers', numWorkers, ...
       'overwriteExisting', overwriteExisting, ...
       'saveData', saveData);

%% Inspect full cross-validation results
[cvResults, bestModels] = getCrossValResults_fa(runIdx, 'baseDir', baseDir);

% Plot cross-validated performance vs estimated dimensionality
plotPerfvsDim_fa(cvResults, ...
                 'bestModels', bestModels);
               
% Find a conservative estimate of optimal total dimensionality for each 
% group.
cutoffPC = 0.95;
numGroups = length(yDims);
xDim_total_fa = nan(1,numGroups);
for groupIdx = 1:numGroups
    xDim_total_fa(groupIdx) ...
        = findSharedDimCutoff_fa(cvResults{groupIdx}(bestModels(groupIdx)).estParams, ...
                                 cutoffPC, ...
                                 'plotSpec', true);
end

%% ================================================================
% 2c) Use the FA estimates above as upper bounds on the search range for
%     optimal DLAG dimensionalities.
%  ======================================================================

% Change other input arguments as appropriate
runIdx = 2;
numFolds = 4;
maxIters = 100; % Limit EM iterations during cross-validation for speedup
% Largest across-group dimensionality we'll consider
xDim_across_max = min(xDim_total_fa);
% All across-group dimensionalities we'll consider
xDims_across = 0:xDim_across_max; 
% Use within-group dimensionalities such that 
% xDim_across + xDim_within = xDim_total_fa
xDims_grid = nan(length(xDims_across),numGroups+1);
for modelIdx = 1:length(xDims_across)
    xDim_across = xDims_across(modelIdx);
    xDim_within = xDim_total_fa - xDim_across;
    xDims_grid(modelIdx,:) = [xDim_across xDim_within];
end

fit_dlag(runIdx, Ytrain, ...
         'baseDir', baseDir, ...
         'method', method, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'xDims_grid', xDims_grid, ...
         'yDims', yDims, ...
         'rGroups', rGroups,...
         'startTau', startTau, ...
         'segLength', segLength, ...
         'init_method', init_method, ...
         'learnDelays', learnDelays, ...
         'maxIters', maxIters, ...
         'freqLL', freqLL, ...
         'freqParam', freqParam, ...
         'minVarFrac', minVarFrac, ...
         'parallelize', parallelize, ...
         'randomSeed', randomSeed, ...
         'numWorkers', numWorkers, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);
     
%% Inspect cross-validation results
% Retrieve cross-validated results for all models in the results directory
[cvResults, ~] = getCrossValResults_dlag(runIdx, 'baseDir', baseDir);

% Plot a variety of performance metrics among the candidate models.
plotPerfvsDim_dlag(cvResults, 'xDims_grid', xDims_grid);

% Select the model with the optimal number of across-group latents
% NOTE: If you limited the number of EM iterations during cross-validation,
%       then it's a good idea to re-fit the optimal model to full
%       convergence. We won't do that here, for concision.
bestModel = getNumAcrossDim_dlag(cvResults, xDims_grid, 'metric', 'R2', 'joint', true);

%% ==========================================================
% 3a) Analyze the spectra of the DLAG loading matrices
%  ================================================================

% Inspect the performance behavior of the orthonormalized DLAG model.
plotOrthPerfvsDim_dlag(bestModel);
dim = getNumOrthDim_dlag(bestModel, 'metric', 'R2')

% Determine the minimum number of latent dimensions needed
% to explain a certain proportion of shared variance, in the joint 
% population, across groups, and for each group individually.
cutoffPC = 0.95;
d_shared = findSharedDimCutoff_dlag(bestModel.estParams, ...
                                    cutoffPC, ...
                                    'plotSpec', true)

%% ==========================================================
% 3b) Evaluate the prominence of individual within- and across-
%     group latent variables.
%  ================================================================

% Save all bootstrap results to a file
boot_fname = sprintf('%s/mat_results/run%03d/bootstrapResults', baseDir, runIdx);

numBootstrap = 100;  % Number of bootstrap samples
alpha = 0.05;        % Construct (1-alpha) bootstrapped confidence intervals 

% Prominence evaluated for all groups jointly
dimGroups_across = num2cell(1:bestModel.xDim_across);
dimGroups_within = {num2cell(1:bestModel.xDim_within(1)), ...
                    num2cell(1:bestModel.xDim_within(2))};
prom = bootstrapProminence(Ytrain, ...
                           bestModel.estParams, ...
                           dimGroups_across, ...
                           dimGroups_within, ...
                           numBootstrap, ...
                           'alpha', alpha, ...
                           'parallelize', parallelize, ...
                           'numWorkers', numWorkers);
save(boot_fname, 'prom');
plotProminence(prom, 'metric', 'LL');

% Prominence evaluated for each group individually
groupProm = groupBootstrapProminence(Ytrain, ...
                                     bestModel.estParams, ...
                                     dimGroups_across, ...
                                     dimGroups_within, ...
                                     numBootstrap, ...
                                     'alpha', alpha, ...
                                     'parallelize', parallelize, ...
                                     'numWorkers', numWorkers);
save(boot_fname, 'groupProm', '-append');
plotGroupProminence(groupProm, 'metric', 'LL');

%% ==========================================================
% 3c) Evaluate how signficantly each set of across-group delays
%     deviates from zero.
%  ================================================================

sig = bootstrapDelaySignificance(Ytrain, ...
                                bestModel.estParams, ...
                                numBootstrap, ...
                                'alpha', alpha, ...
                                'parallelize', parallelize, ...
                                'numWorkers', numWorkers);
save(boot_fname, 'sig', '-append');
plotDelaySignificance(sig, bestModel.estParams, ...
                      binWidth, bestModel.rGroups);

%% ==========================================================
% 3d) Construct bootstrapped confidence intervals for latent delays
%     and timescales.
%  ================================================================

bootParams = bootstrapGPparams(Ytrain, ...
                               bestModel.estParams, ...
                               binWidth, ...
                               numBootstrap, ...
                               'alpha', alpha, ...
                               'parallelize', parallelize, ...
                               'numWorkers', numWorkers, ...
                               'segLength', Inf, ...
                               'tolLL', 1e-4);
save(boot_fname, 'bootParams', '-append');
plotBootstrapGPparams_dlag(bestModel.estParams, bootParams, ...
                           binWidth, bestModel.rGroups, 'overlayParams', true);