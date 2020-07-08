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

%% Define DLAG ground truth model parameters

rng('shuffle');

% Dataset size characteristics
N = 200;                            % Total number of sequences (trials)
numTrain = 100;                     % Define train and test splits
numTest = N - numTrain;
T = 25;                             % Number of samples per sequence
binWidth = 20;                      % Sample period of ground truth data (take only one sample every binWidth samples)
yDims = [10 10];                    % Dimensionalities of each observed group
numGroups = length(yDims);          % Total number of groups
% Used for converting full C matrix into Cs cell array
block_idxs = get_block_idxs(yDims);
xDim_across = 4;                    % Across-group latent dimensionality
xDim_within = [2 2];                % Within-group latent dimensionalities
snr = [0.5 0.5];                    % Signal-to-noise ratios of each group
centered = false;                   % Let d be non-zero

% GP Timescales
min_tau = 30;                       % Lower-bound of GP timescale range
max_tau = 130;                      % Upper-bound of GP timescale range

% GP noise variances                
min_eps = 1e-5;                     % Lower-bound of GP noise range
max_eps = 1e-5;                     % Upper-bound of GP noise range

% Delays
min_delay = -25;                 % Lower-bound of delay range, in samples (time steps)
max_delay = 25;                  % Uppder-bound of delay range, in samples (time steps)

%% Randomly generate data from a DLAG model
[Ys, Xs, params] = simdata_dlag(N, T, binWidth, yDims, ...
                                xDim_across, xDim_within, ...
                                min_tau, max_tau, min_eps, max_eps, ...
                                min_delay, max_delay, snr, centered);

% Create train and test splits.
% For observations, first convert the data into a format compatible with 
% fit_dlag.m
Y = cat(1, Ys{:});
Y = dat2seq(Y);
Ytrain = Y(1:numTrain);
Ytest = Y(numTrain+1:end);

% Do the same for ground truth latent trajectories, so data is compatible
% with DLAG plotting functions
X = cat(1, Xs{:});
X = dat2seq(X,'datafield','xgt');
Xtrain = X(1:numTrain);
Xtest = X(numTrain+1:end);

%% Save some of the ground truth parameters that we care about
trueParams.C = blkdiag(params.Cs{:});
trueParams.R = blkdiag(params.Rs{:});
trueParams.d = cat(1, params.ds{:});
trueParams.gamma_across = (binWidth./(params.taus_across)).^2;
trueParams.eps_across = params.eps_across;
for groupIdx = 1:numGroups
    trueParams.gamma_within{groupIdx} = (binWidth./params.taus_within{groupIdx}).^2;
    trueParams.eps_within{groupIdx} = params.eps_within{groupIdx};
end
trueParams.DelayMatrix = cat(1, params.delays{:})./binWidth;
trueParams.covType = 'rbf';
trueParams.notes.RforceDiagonal = true;
trueParams.notes.learnKernelParams = true;
trueParams.notes.learnGPNoise = false;
trueParams.xDim_across = xDim_across;
trueParams.xDim_within = xDim_within;
trueParams.yDims = yDims;

%% Save generated data, along with ground truth parameters
save('mat_sample/dlag_demo_data_synthetic.mat', 'trueParams', 'Ys', 'Ytrain', 'Ytest', 'Xs',...
     'Xtrain', 'Xtest', 'snr', 'binWidth');