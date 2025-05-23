% =================================================================
% DLAG-frequency DEMO: Run this script from the main DLAG directory
% =================================================================
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% ================
% 0) Load demo data 
% =================

% Synthetic data generated from a DLAG model
dat_file = 'demo/data/dlag-freq_demo_data';
fprintf('Reading from %s \n',dat_file);
load(dat_file);

% Color scheme
GTCOLOR = '#7F7F7F';
TIMECOLOR = 'k';
FREQCOLOR = '#D35FBC';

%% ===================================================================
% 1a) Compare the runtimes of time domain and frequency domain fitting
% ====================================================================

% Common arguments across both approaches

% Let's explicitly define all of the optional arguments, for 
% the sake of demonstration:
runIdx = 7;               % Results will be saved in baseDir/mat_results/runXXX/,  
                          % where XXX is runIdx. Use a new runIdx for each dataset.
baseDir = './demo';       % Base directory where results will be saved
overwriteExisting = true; % Control whether existing results files are overwritten
saveData = false;         % Set to true to save train and test data (not recommended)
binWidth = 20;            % Sample period / spike count bin width, in units of time (e.g., ms)
numFolds = 0;             % Number of cross-validation folds (0 means no cross-validation)
xDims_across = 2;         % This number of across-group latents matches the synthetic ground truth
xDims_within = {0, 0};    % These numbers match the within-group latents in the synthetic ground truth
yDims = [10 10];          % Number of observed features (neurons) in each group (area)
covType = 'rbf';          % Type of GP covariance kernel ('rbf' or 'sg')
rGroups = [1 2];          % For performance evaluation, we can regress group 2's activity with group 1
startTau = 2*binWidth;    % Initial timescale, in the same units of time as binWidth
segLength = Inf;          % Largest trial segment length, in no. of time points
init_method = 'static';   % Initialize DLAG with fitted pCCA parameters
learnDelays = true;       % Set to false if you want to fix delays at their initial value
maxIters = 10;            % Limit the number of EM iterations (not recommended for final fitting stage)
tolLL = 1e-8;             % Log-likelihood convergence tolerance
freqLL = 1;               % Check for data log-likelihood convergence every freqLL EM iterations
freqParam = 10;           % Store intermediate delay and timescale estimates every freqParam EM iterations
minVarFrac = 0.001;       % Private noise variances will not be allowed to go below this value
verbose = true;           % Toggle printed progress updates
randomSeed = 0;           % Seed the random number generator, for reproducibility

%% =====================================================================
% 1b) Time domain fitting
%     NOTE: You can skip this section if a trained model already exists.
% ======================================================================

method = 'dlag';  % Time domain 'dlag' fitting

fit_dlag(runIdx, seqTrue, ...
         'baseDir', baseDir, ...
         'method', method, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'xDims_across', xDims_across, ...
         'xDims_within', xDims_within, ...
         'yDims', yDims, ...
         'covType', covType, ...
         'rGroups', rGroups,...
         'startTau', startTau, ...
         'segLength', segLength, ...
         'init_method', init_method, ...
         'learnDelays', learnDelays, ...
         'maxIters', maxIters, ...
         'tolLL', tolLL, ...
         'freqLL', freqLL, ...
         'freqParam', freqParam, ...
         'minVarFrac', minVarFrac, ...
         'parallelize', false, ...
         'verbose', verbose, ...
         'randomSeed', randomSeed, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);

%% =====================================================================
% 1c) Frequency domain fitting
%     NOTE: You can skip this section if a trained model already exists.
% ======================================================================

method = 'dlag-freq';  % Frequency domain 'dlag-freq' fitting

fit_dlag(runIdx, seqTrue, ...
         'baseDir', baseDir, ...
         'method', method, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'xDims_across', xDims_across, ...
         'xDims_within', xDims_within, ...
         'yDims', yDims, ...
         'covType', covType, ...
         'rGroups', rGroups,...
         'startTau', startTau, ...
         'segLength', segLength, ...
         'init_method', init_method, ...
         'learnDelays', learnDelays, ...
         'maxIters', maxIters, ...
         'tolLL', tolLL, ...
         'freqLL', freqLL, ...
         'freqParam', freqParam, ...
         'minVarFrac', minVarFrac, ...
         'parallelize', false, ...
         'verbose', verbose, ...
         'randomSeed', randomSeed, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);

%% =======================
% 1d) Investigate runtimes
% ========================

xDim_across = xDims_across;
xDim_within = cell2mat(xDims_within);
fprintf('Avg. runtime per iteration (s):\n')

% Time domain
method = 'dlag';
res = getModel_dlag(runIdx, xDim_across, xDim_within, ...
                    'baseDir', baseDir, 'method', method);
fprintf('    time domain:    %f\n', mean(res.iterTime))

% Frequency domain
method = 'dlag-freq';
res = getModel_dlag(runIdx, xDim_across, xDim_within, ...
                    'baseDir', baseDir, 'method', method);
fprintf('    freq domain:    %f\n', mean(res.iterTime))

%% ======================================================================
% 2a) Fully fit a DLAG model using the frequency domain approach
%     NOTE: You can skip to Section 2b if a trained model already exists.
% =======================================================================

method = 'dlag-freq';     % Time domain 'dlag' or frequency domain 'dlag-freq' fitting

runIdx = 8;               % Results will be saved in baseDir/mat_results/runXXX/,  
                          % where XXX is runIdx. Use a new runIdx for each dataset.
baseDir = './demo';       % Base directory where results will be saved
overwriteExisting = true; % Control whether existing results files are overwritten
saveData = false;         % Set to true to save train and test data (not recommended)
binWidth = 20;            % Sample period / spike count bin width, in units of time (e.g., ms)
numFolds = 0;             % Number of cross-validation folds (0 means no cross-validation)
xDims_across = 2;         % This number of across-group latents matches the synthetic ground truth
xDims_within = {0, 0};    % These numbers match the within-group latents in the synthetic ground truth
yDims = [10 10];          % Number of observed features (neurons) in each group (area)
covType = 'rbf';          % Type of GP covariance kernel ('rbf' or 'sg')
rGroups = [1 2];          % For performance evaluation, we can regress group 2's activity with group 1
startTau = 2*binWidth;    % Initial timescale, in the same units of time as binWidth
segLength = Inf;          % Largest trial segment length, in no. of time points
init_method = 'static';   % Initialize DLAG with fitted pCCA parameters
learnDelays = true;       % Set to false if you want to fix delays at their initial value
maxIters = 1e4;           % Limit the number of EM iterations (not recommended for final fitting stage)
tolLL = 1e-8;             % Log-likelihood convergence tolerance
freqLL = 1;               % Check for data log-likelihood convergence every freqLL EM iterations
freqParam = 10;           % Store intermediate delay and timescale estimates every freqParam EM iterations
minVarFrac = 0.001;       % Private noise variances will not be allowed to go below this value
verbose = true;           % Toggle printed progress updates
randomSeed = 0;           % Seed the random number generator, for reproducibility

fit_dlag(runIdx, seqTrue, ...
         'baseDir', baseDir, ...
         'method', method, ...
         'binWidth', binWidth, ...
         'numFolds', numFolds, ...
         'xDims_across', xDims_across, ...
         'xDims_within', xDims_within, ...
         'yDims', yDims, ...
         'covType', covType, ...
         'rGroups', rGroups,...
         'startTau', startTau, ...
         'segLength', segLength, ...
         'init_method', init_method, ...
         'learnDelays', learnDelays, ...
         'maxIters', maxIters, ...
         'tolLL', tolLL, ...
         'freqLL', freqLL, ...
         'freqParam', freqParam, ...
         'minVarFrac', minVarFrac, ...
         'parallelize', false, ...
         'verbose', verbose, ...
         'randomSeed', randomSeed, ...
         'overwriteExisting', overwriteExisting, ...
         'saveData', saveData);

%% ========================
% 2b) Check fitting results
% =========================

% Retrieve the fitted model of interest
xDim_across = xDims_across;
xDim_within = cell2mat(xDims_within);
res = getModel_dlag(runIdx, xDim_across, xDim_within, ...
                    'baseDir', baseDir, 'method', method);

% Plot training progress of various quantities
plotFittingProgress(res, ...
                    'freqLL', freqLL, ...
                    'freqParam', freqParam, ...
                    'units', 'ms');

%% ======================================
% 2c) Visualize recovery of GP parameters
% =======================================

% Plot estimated and ground truth delays and across-group GP timescales
% together on the same plot.
plotGPparams_withGT_dlag(res.estParams, trueParams, res.binWidth,...
                         res.rGroups, 'units', 'ms');

%% ============================================
% 2d) Visualize recovery of latent time courses
% =============================================

% NOTE: We fit a DLAG model via the frequency domain. Here, we'll use
%       those fitted parameters, but perform all inference in the time
%       domain.

trialIdx = 1;  % Choose an example trial to plot

% Estimated latents
% NOTE: In general, latent time courses are unordered, and will not match
%       the order in seqTrue.
[seqEst, ~] = exactInferenceWithLL_dlag(seqTrue, res.estParams);
plotDimsVsTime_dlag(seqEst(trialIdx), 'xsm', res.estParams, res.binWidth, ...
                  'nPlotMax', 1, ...
                  'plotSingle', false, ...
                  'plotMean', true, ...
                  'units', [], ...
                  'trialGroups', {[1]}, ...
                  'trialColors', {FREQCOLOR});

% Ground truth latents
plotDimsVsTime_dlag(seqTrue(trialIdx), 'xsm', trueParams, res.binWidth, ...
                  'nPlotMax', 1, ...
                  'plotSingle', false, ...
                  'plotMean', true, ...
                  'units', [], ...
                  'trialGroups', {[1]}, ...
                  'trialColors', {GTCOLOR});

%% =====================================================================
% 3) Using the same set of parameters, compare time and frequency domain
%    inference
% ======================================================================

trialIdx = 1; % Choose an example trial to plot

% Load fitted parameters from Section 2
runIdx = 8;
res = getModel_dlag(runIdx, xDim_across, xDim_within, ...
                    'baseDir', baseDir, 'method', method);

fprintf('\nLatent inference runtime (s):\n');

% Time domain estimate
tic;
[seqEst_time, ~] = exactInferenceWithLL_dlag(seqTrue, res.estParams);
runtime = toc;
fprintf('    time domain:  %1.4f\n', runtime);

% Frequency domain estimate
% Take the FFT of each neuron's activity on each trial
seqTrue = fftseq(seqTrue,'y','yfft');
% Perform inference in the frequency domain
tic;
[seqEst_freq, ~, ~] = exactInferenceWithLL_freq(seqTrue, res.estParams);
runtime = toc;
fprintf('    freq domain:  %f\n', runtime);
% Convert estimated latents back to the time domain
seqEst_freq = freq2time_dlag(seqEst_freq, res.estParams, ...
    'infield', 'xfft', 'outfield', 'xsm');

% Let's plot ground truth, time and frequency domain estimates on top of
% each other
latentIdx = 1;      % Choose one latent to plot
reorder = [1 2];    % Reorder estimated latents as needed
rescale = [-1 -1];   % Flip estimated latents as needed

% Notice the edge effects that appear in the frequency domain estimates
figure;
hold on;
% Ground truth
h1 = plot(1:length(seqTrue(trialIdx).xsm(latentIdx,:)), ...
          seqTrue(trialIdx).xsm(latentIdx,:), ...
          'color', GTCOLOR', ...
          'linestyle', '-', ...
          'linewidth', 2.0);
% Time domain approach
h2 = plot(1:length(seqEst_time(trialIdx).xsm(reorder(latentIdx),:)), ...
          seqEst_time(trialIdx).xsm(reorder(latentIdx),:)*rescale(latentIdx), ...
          'color', TIMECOLOR', ...
          'linestyle', '-', ...
          'linewidth', 2.0);
% Frequency domain approach
h3 = plot(1:length(seqEst_freq(trialIdx).xsm(reorder(latentIdx),:)), ...
          seqEst_freq(trialIdx).xsm(reorder(latentIdx),:)*rescale(latentIdx), ...
          'color', FREQCOLOR', ...
          'linestyle', '-', ...
          'linewidth', 2.0);
legend([h1, h2, h3], {'Ground truth', 'Time', 'Freq'});
xlabel('Time point');
ylabel('x');

% Error between time and frequency domain approaches as a function of time.
% Again, note the error increase at the edges.
Xtime = cat(3,seqEst_time.xsm);
Xfreq = cat(3,seqEst_freq.xsm);
Xerr = sqrt(mean((Xtime(reorder,:,:) - Xfreq(reorder,:,:)).^2,3)); % RMSE
figure;
for j = 1:size(Xerr,1)
    subplot(size(Xerr,1),1,j);
    hold on;
    plot(Xerr(j,:),'k-', 'linewidth', 1.5);
    xlabel('Time');
    ylabel('Error');
    title(sprintf('Latent %d', j));
end
