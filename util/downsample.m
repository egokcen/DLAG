function Yout = downsample(Yin, period)
%
% Yout = downsample(Yin, period)
%
% Description: Drop samples from a set of sequences, Yin, with a new sample
%              rate of period.
%
% Arguments:
%     Yin    -- (1 x numGroups) cell array; list of data matrices 
%               {(y1Dim x T x N), (y2Dim x T x N), ...}
%     period -- int; downsample period, in units of samples
%
% Outputs:
%     Yout   -- (1 x numGroups) cell array; list of downsampled data matrices 
%               {(y1Dim x Tdown x N), (y2Dim x Tdown x N), ...}
%
% Author: Evren Gokcen

    numGroups = length(Yin);
    Yout = cell(size(Yin));
    T = size(Yin{1},2);
    for groupIdx = 1:numGroups
        sampleIdxs = 1:period:T;
        Yout{groupIdx} = Yin{groupIdx}(:,sampleIdxs,:);
    end
end