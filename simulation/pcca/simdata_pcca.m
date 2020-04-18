function [Ys, X, params] = simdata_pcca(N, yDims, xDim, snr, centered)
% 
% [Ys, X, params] = simdata_pcca(N, yDims, xDim, snr, centered) 
%
% Description: Generate simulated data according to a pCCA model.
%
% Arguments:
%     N        -- int; number of data points to generate
%     yDims    -- (1 x numGroups) array; List of dimensionalities of
%                 observed data, [y1Dim, y2Dim, ...]
%     xDim     -- int; Dimensionality of latents, X
%     snr      -- (1 x numGroups) array; List of signal-to-noise ratios,
%                 defined as trace(CC') / trace(R)
%     centered -- boolean; true if data is to be zero-centered
%
% Outputs:
%     Ys        -- (1 x numGroups) cell array; list of data matrices 
%                  {(y1Dim x N), (y2Dim x N), ...}
%     X         -- (xDim x N) array; latent data
%     params.Cs -- (1 x numGroups) cell array; List of factor loadings 
%                  {(y1Dim x xDim), (y2Dim x xDim), ...}
%     params.Rs -- (1 x numGroups) cell array; Blocks of uniqueness 
%                  matrix {(y1Dim x y1Dim), (y2Dim x y2Dim), ...}
%     params.ds -- (1 x numGroups) cell array; List of data means 
%                  {(y1Dim x 1), (y2Dim x 1), ...}
%
% Author: Evren Gokcen

    numGroups = length(yDims);
    
    % Generate latent data
    X = randn(xDim, N);
    
    % Generate PCCA model parameters with specified SNR
    Cs = cell(1, numGroups);
    Rs = cell(1, numGroups);
    ds = cell(1, numGroups);
    Ys = cell(1, numGroups);
    for groupIdx = 1:numGroups
        yDim = yDims(groupIdx);
        % Loading matrices
        C = randn(yDim, xDim);
        % Uniqueness matrices
        R = randn(yDim, yDim);
        % Ensure that each R is a valid covariance matrix
        R = R * R';
        
        % Enforce desired signal-to-noise ratios
        varCC = trace(C * C');
        varNoise_desired = varCC / snr(groupIdx);
        varNoise_current = trace(R);
        R = R .* (varNoise_desired / varNoise_current);
        
        % Center the data, if desired
        if centered
            d = zeros(yDim,1);
        else
            d = randn(yDim,1);
        end
        
        % Generate observations
        ns = mvnrnd(zeros(1,yDim), R, N)';
        Y = C * X + repmat(d, 1, N) + ns;
        
        % Collect all outputs
        Cs{groupIdx} = C;
        Rs{groupIdx} = R;
        ds{groupIdx} = d;
        Ys{groupIdx} = Y;
        
    end
    
    params.Cs = Cs;
    params.Rs = Rs;
    params.ds = ds;

end