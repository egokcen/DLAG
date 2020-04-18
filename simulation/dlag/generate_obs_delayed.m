function [Ys, X_clipped] = generate_obs_delayed(X, Cs, Rs, ds, max_delay, delays)
%
% [Ys, X_clipped] = generate_obs_delayed(X, Cs, Rs, ds, max_delay, delays)
%
% Description: Given a sequence of latent variables, X, generate observed
%              data according to a linear observation model.
%
% Arguments:
%     X  -- (xDim x T) array; simulated latent sequences
%     Cs -- (1 x numGroups) cell array; List of factor loadings 
%           {(y1Dim x xDim), (y2Dim x xDim), ...}
%     Rs -- (1 x numGroups) cell array; Blocks of uniqueness 
%           matrix {(y1Dim x y1Dim), (y2Dim x y2Dim), ...}
%     ds -- (1 x numGroups) cell array; List of data means 
%           {(y1Dim x 1), (y2Dim x 1), ...}
%     max_delay -- int; maximum delay (+ or -), in ms
%     delays -- (1 x numGroups) cell array; Delays between
%               latents and each group. delays{i} (1 x xDim) array.
%     
% Outputs:
%     Ys -- (1 x numGroups) cell array; list of data matrices 
%           {(y1Dim x T), (y2Dim x T), ...}
%     X_clipped -- (1 x numGroups) cell array; list of latent data
%                  matrices, clipped to be the same length as corresponding
%                  observations in Ys
%                  {(x1Dim x T), (x2Dim x T), ...}
%
% Author: Evren Gokcen

    [xDim, T] = size(X);
    Tnew = T - 2*max_delay;
    numGroups = length(Cs);
    Ys = cell(1, numGroups);
    X_clipped = cell(1, numGroups);
    for groupIdx = 1:numGroups
        yDim = size(Rs{groupIdx},1);
        % Delay each latent separately, relative to the current group
        X_delay = zeros(xDim, Tnew);
        curr_delays = delays{groupIdx};
        for l = 1:xDim
            timeIdxs = time_shift(T, curr_delays(l), max_delay);
            X_delay(l,:) = X(l,timeIdxs);
        end
        ns = mvnrnd(zeros(1,yDim), Rs{groupIdx}, Tnew)';
        Ys{groupIdx} = Cs{groupIdx} * X_delay + repmat(ds{groupIdx}, 1, Tnew) + ns;
        X_clipped{groupIdx} = X_delay;
    end
end