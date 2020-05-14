function seq = dat2seq(dat, varargin)
%
% seq = dat2seq(dat,...) 
%
% Description: Convert data matrix of sequences and trials into format
%              compatible with DLAG functions.
%
% Arguments:
%
%     Required:
%
%     dat -- (D x T x N) array; N sequences of length T and dimensionality 
%            D
%
%     Optional:
%     
%     datafield -- string; Name of data field in seq (default: 'y')
%
% Outputs:
%     seq -- structure whose nth entry (corresponding to the nth sequence)
%            has fields
%                    trialId   -- unique trial identifier  
%                    T (1 x 1) -- number of timesteps
%                    y (D x T) -- continuous valued data 
%
% Author: Evren Gokcen

    datafield = 'y';
    extraOpts = assignopts(who, varargin);
    [D, T, N] = size(dat);
    for n = 1:N
        seq(n).trialId = n;
        seq(n).T = T;
        seq(n).(datafield) = dat(:,:,n);
    end
end