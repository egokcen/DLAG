% dlag_groundtruth_data.m
%
% Description: This script generates synthetic datasets from the DLAG
%              model. Ground truth data and parameters are saved to a file.
%              Data is formatted to be compatible with fit_dlag.m
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu

%% Define DLAG ground truth model parameters

rng('shuffle');

% Dataset size characteristics
N = 100;                          % Number of sequences (trials)
T = 25;                           % Number of samples per sequence
binWidth = 20;                    % Sample period of ground truth data
yDims = [10 10];                  % Dimensionalities of each observed group
numGroups = length(yDims);        % Total number of groups
xDim_across = 4;                  % Across-group latent dimensionality
xDim_within = [2 2];              % Within-group latent dimensionalities
snr = [1.0 1.0];                  % Signal-to-noise ratios of each group
covType = 'rbf';                  % Specify GP covariance type

% GP Timescales
param_lims.tau_lim = [20 150];     % GP timescale range

% GP noise variances                
param_lims.eps_lim = [1e-5 1e-5];  % GP noise range

% Delays
param_lims.delay_lim = [-30 30];   % Delay range, in samples

% Center frequencies
param_lims.nu_lim = [0 10]./1000;  % Center frequency range, in 1/ms (to
                                   % match units of delays and timescales)

%% Randomly generate data from a DLAG model

[seqTrue, trueParams] = simdata_dlag(N, T, binWidth, yDims, ...
                                     xDim_across, xDim_within, snr, ...
                                     covType, param_lims);   

% Add relevant notes to trueParams                        
trueParams.notes.RforceDiagonal = true;
trueParams.notes.learnKernelParams = true;
trueParams.notes.learnGPNoise = false;

%% Visualize the ground truth

% GP parameters
plotGPparams_dlag(trueParams, binWidth, [1 2], ...
                  'plotAcross', true, ...
                  'plotWithin', true, ...
                  'units', 'ms');

% Latent timecourses
plotDimsVsTime_dlag(seqTrue, 'xsm', trueParams, binWidth, ...
                   'nPlotMax', 1, ...
                   'plotSingle', false, ...
                   'plotMean', true, ...
                   'units', 'ms');
                    
%% Save generated data, along with ground truth parameters

save('demo/data/dlag-freq_demo_data.mat', 'seqTrue', 'trueParams', ...
     'snr', 'binWidth');