% pcca_groundtruth_data.m
%
% Description: This script generates synthetic datasets from the pCCA
%              model. Ground truth data and parameters are saved to a file.
%              Data is formatted to be compatible with fit_static.m
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     07 Apr 2020 -- Initial full revision.
%     14 May 2020 -- Updated to remain compatible with 'dat' argument
%                    format in fit_dlag.m

%% Define pCCA ground truth model parameters

rng('shuffle');

% Dataset size characteristics
N = 200;                            % Total number of sequences (trials)
numTrain = 100;                     % Define train and test splits
numTest = N - numTrain;
T = 20;                             % Number of samples per sequence
binWidth = 20;                      % Sample period of ground truth data (take only one sample every binWidth samples)
yDims = [10 10];                    % Dimensionalities of each observed group
yDim = sum(yDims);
numGroups = length(yDims);          % Total number of groups
% Used for converting full C matrix into Cs cell array
block_idxs = get_block_idxs(yDims);
xDim = 4;                           % Across-group latent dimensionality
snr = [1.0 1.0];                    % Signal-to-noise ratios of each group
centered = false;                   % Let d be non-zero

%% Randomly generate data from a pCCA model
% pCCA treats all data points the same, regardless of trial or location
% in time
numSamples = N*T; 
[Ys, X, params] = simdata_pcca(numSamples, yDims, xDim, snr, centered);

% Create train and test splits.
% First convert the data into a format compatible with fit_static.m
Y = zeros(yDim,T,N);
Tlist = repmat(T,1,N);
time_blocks = get_block_idxs(Tlist);
group_blocks = get_block_idxs(yDims);
for n = 1:N
    currTrial = time_blocks{n}(1):time_blocks{n}(2);
    for groupIdx = 1:numGroups
        currGroup = group_blocks{groupIdx}(1):group_blocks{groupIdx}(2);
        Y(currGroup,:,n) = Ys{groupIdx}(:,currTrial);
    end
end
Y = dat2seq(Y);
Ytrain = Y(1:numTrain);
Ytest = Y(numTrain+1:end);

% Do the same for ground truth latent trajectories, so data is compatible
% with plotting functions
Xgt = zeros(xDim,T,N);
for n = 1:N
    currTrial = time_blocks{n}(1):time_blocks{n}(2);
    Xgt(:,:,n) = X(:,currTrial);
end
Xgt = dat2seq(Xgt,'datafield','xgt');
Xtrain = Xgt(1:numTrain);
Xtest = Xgt(numTrain+1:end);

%% Save some of the ground truth parameters that we care about
trueParams.Cs = params.Cs;
trueParams.Rs = params.Rs;
trueParams.d = params.ds;
trueParams.xDim = xDim;
trueParams.yDims = yDims;

%% Save generated data, along with ground truth parameters
save('mat_sample/pcca_demo_data_synthetic.mat', 'trueParams', 'Ys',...
     'Ytrain', 'Ytest', 'X', 'Xtrain', 'Xtest', 'snr', 'binWidth');