% dlag_groundtruth_data.m
%
% Description: This script generates synthetic datasets from the DLAG
%              model. Ground truth data and parameters are saved to a file.
%              Data is formatted to be compatible with fit_dlag.m
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     07 Apr 2020 -- Initial full revision.
%     14 May 2020 -- Updated to remain compatible with 'dat' argument
%                    format in fit_dlag.m
%     14 Oct 2020 -- Updated to reflect changes to simdata_dlag.m
%     18 Mar 2021 -- Updated to reflect changes to simdata_dlag.m

%% Define DLAG ground truth model parameters

rng('shuffle');

% Dataset size characteristics
N = 100;                            % Number of sequences (trials)
T = 25;                             % Number of samples per sequence
binWidth = 20;                      % Sample period of ground truth data
yDims = [10 10];                    % Dimensionalities of each observed group
numGroups = length(yDims);          % Total number of groups
xDim_across = 4;                    % Across-group latent dimensionality
xDim_within = [2 2];                % Within-group latent dimensionalities
snr = [0.5 0.5];                    % Signal-to-noise ratios of each group

% GP Timescales
tau_lim = [30 130];                 % GP timescale range

% GP noise variances                
eps_lim = [1e-5 1e-5];              % GP noise range

% Delays
delay_lim = [-25 25];               % Delay range, in samples

%% Randomly generate data from a DLAG model

[seqTrue, trueParams] = simdata_dlag(N, T, binWidth, yDims, ...
                                xDim_across, xDim_within, snr, ...
                                tau_lim, eps_lim, delay_lim);                        
% Add relevant notes to trueParams                        
trueParams.notes.RforceDiagonal = true;
trueParams.notes.learnKernelParams = true;
trueParams.notes.learnGPNoise = false;

% Investigate generated GP parameters
plotGPparams_dlag(trueParams, binWidth, [1 2], ...
                  'plotAcross', true, ...
                  'plotWithin', true, ...
                  'units', 'ms');

%% Save generated data, along with ground truth parameters
save('mat_sample/dlag_demo_data_synthetic.mat', 'seqTrue', 'trueParams', ...
     'snr', 'binWidth');