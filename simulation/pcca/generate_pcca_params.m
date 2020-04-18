function [Cs, Rs, ds] = generate_pcca_params(yDims, xDim, snr, centered)
% 
% [Cs, Rs, ds] = generate_pcca_params(yDims, xDim, snr, centered)
%
% Description: Generate parameters that define a pCCA model with a 
%              specified signal-to-noise ratios.
% 
% Arguments:
%     yDims    -- (1 x numGroups) array; List of dimensionalities of
%                 observed data, [Y1, Y2,...]
%     xDim     -- int; Dimensionality of latents, X
%     snr      -- (1 x numGroups) array; List of signal-to-noise ratios,
%                 defined as trace(CC') / trace(R)
%     centered -- boolean; True if data is to be zero-centered
%
% Outputs:
%     Cs -- (1 x numGroups) cell array; List of factor loadings {(y1Dim x
%           xDim), (y2Dim x xDim), ...}
%     Rs -- (1 x numGroups) cell array; Blocks of uniqueness matrix {(y1Dim
%           x y1Dim), (y2Dim x y2Dim), ...}
%     ds -- (1 x numGroups) cell array; List of data means {(y1Dim x 1),
%           (y2Dim x 1), ...}
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     07 Apr 2020 -- Initial full revision.

    numGroups = length(yDims);
    % Initialize output arrays
    Cs = cell(1, numGroups);
    Rs = cell(1, numGroups);
    ds = cell(1, numGroups);
    for groupIdx = 1:numGroups
        C = randn(yDims(groupIdx), xDim);
        R = randn(yDims(groupIdx), yDims(groupIdx));
        % Ensure that R is a valid covariance matrix
        R = R * R';
        
        % Enforce the desired signal-to-noise ratios
        varCC = trace(C * C');
        varNoise_desired = varCC / snr(groupIdx);
        varNoise_current = trace(R);
        R = R .* (varNoise_desired / varNoise_current);
        
        % Center the data, if desired
        if centered
            d = zeros(yDims(groupIdx),1);
        else
            d = randn(yDims(groupIdx),1);
        end
        
        % Collect all outputs
        Cs{groupIdx} = C;
        Rs{groupIdx} = R;
        ds{groupIdx} = d;
        
    end

end