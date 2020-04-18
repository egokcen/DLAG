function precomp = makePrecomp_timescales(seq,params)
%
% [precomp] = makePrecomp_timescales(seq,params)
%
% Description: Precompute posterior covariances for the DLAG fitting 
%              procedure.
%              NOTE: It might be a good idea to provide a mex 
%                    implementation of this in the future.
%
% Arguments:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                     trialId      -- unique trial identifier
%                     T (1 x 1)    -- number of timesteps
%                     y (yDim x T) -- neural data
%
%     params -- Structure containing DLAG ACROSS-GROUP model parameters. 
%               Contains the fields
% 
%               covType -- string; type of GP covariance (e.g., 'rbf')
%               gamma   -- (1 x xDim) array; GP timescales
%                          in ms are given by 'stepSize ./ sqrt(gamma)'                                                    
%               eps     -- (1 x xDim) GP noise variances
%               d            -- (yDim x 1) array; observation mean
%               C            -- (yDim x (numGroups*xDim)) array;
%                               mapping between low- and high-d spaces
%               R            -- (yDim x yDim) array; observation noise
%                               covariance 
%               DelayMatrix  -- (numGroups x xDim_across) array;
%                               delays from across-group latents to 
%                               observed variables. NOTE: Delays are
%                               reported as (real-valued) number of
%                               time-steps.
%               xDim    -- int; number of across-group latent variables
%               yDims   -- (1 x numGroups) array; 
%                          dimensionalities of each observed group
%
% Outputs
%     precomp -- Structure whose ith entry contains precomputations for 
%                the i-th across-group state:
%                Tdif -- (numGroups*T x numGroups*T) array; matrix of time
%                        differences.
%                Tall -- (1 x N) array; List of trial lengths
%                params -- structure of the same format as 'params' above.
%                Tu   -- structure whose jth entry, corresponding to a
%                        group of trials of the same length, contains 
%                        the following:
%                        nList -- List of trial IDs belonging to this group
%                        T     -- int; Length of all trials in this group
%                        numTrials -- int; Number of trials in this group
%                        PautSum -- (numGroups*T x numGroups*T) array;
%                                   Precomputed posterior GP covariance 
%                                   based on the observations in this group            
% 
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     18 Mar 2020 -- Initial full revision. 

% Initialize relevant variables
xDim = params.xDim;
yDims = params.yDims;
numGroups = length(yDims);
Tall = [seq.T];
Tmax = max(Tall);
Tdif = repmat(1:Tmax,numGroups,1);
Tdif = repmat(Tdif(:)',numGroups*Tmax,1) - repmat(Tdif(:),1,numGroups*Tmax);

% Assign some helpful precomp items
precomp(xDim).Tdif = Tdif;
for i = 1:xDim
    precomp(i).Tdif = Tdif;
    %precomp(i).absDif = abs(Tdif);
    precomp(i).Tall   = Tall;
    precomp(i).params = params;
    
end
% Find unique numbers of trial lengths
Tu = unique(Tall);
% Loop once for each state dimension (each GP)
for i = 1:xDim
    for j = 1:length(Tu)
        T = Tu(j);
        precomp(i).Tu(j).nList = find(Tall == T);
        precomp(i).Tu(j).T = T;
        precomp(i).Tu(j).numTrials = length(precomp(i).Tu(j).nList);
        precomp(i).Tu(j).PautoSUM  = zeros(T*numGroups);
    end
end

% Fill out PautoSum
% Loop once for each state dimension (each GP)
for i = 1:xDim
    % Loop once for each trial length (each of Tu)
    for j = 1:length(Tu)
        % Loop once for each trial (each of nList)
        for n = precomp(i).Tu(j).nList
            xsm_i = seq(n).xsm(i:xDim:xDim*numGroups,:);
            xsm_i = reshape(xsm_i,1,numGroups*Tu(j));
            
            precomp(i).Tu(j).PautoSUM = precomp(i).Tu(j).PautoSUM +...
                seq(n).VsmGP(:,:,i) +...
                xsm_i' * xsm_i;
        end
    end
end
