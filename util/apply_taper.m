function seq_out = apply_taper(seq_in, taper, datafield)
%
% seq_out = apply_taper(seq_in, taper, datafield)
%
% Apply a taper to each trial of seq_in(n).datafield. Importantly,
% rescale the tapered data so that the mean and total variance of the data 
% is conserved.
%
% Arguments:
%
%     seq_in    -- data structure, whose nth entry (corresponding to
%                  the nth trial) has fields
%                      trialId             -- unique trial identifier
%                      T (1 x 1)           -- number of timesteps
%                      datafield (dim x T) -- sequential data
%     taper     -- One of Matlab's builtin window functions that returns a
%                  vector of tapering weights. See, for example, Matlab's 
%                  'hamming' or 'hann' functions.
%     datafield -- string; field name for the data to be tapered
%
% Outputs:
%
%     seq_out   -- data structure, whose nth entry (corresponding to
%                  the nth trial) has fields
%                      trialId             -- unique trial identifier
%                      T (1 x 1)           -- number of timesteps
%                      datafield (dim x T) -- tapered sequential data
%
% Examples:
%     >> seq_out = apply_taper(seq_in, @(T) hamming(T, 'periodic'), 'y');
%     >> seq_out = apply_taper(seq_in, @(T) tukeywin(T, 0.1), 'y');
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     14 Aug 2024 -- Initial full revision.

% Get the mean and variance of the raw data
var_raw = var([seq_in.(datafield)],0,2);  % Raw variance
mean_raw = mean([seq_in.(datafield)],2);  % Raw mean

seq_out = [];
for n = 1:length(seq_in)
    seq_out(n).trialId = seq_in(n).trialId;
    seq_out(n).T = seq_in(n).T;
    % Apply the periodic taper to normalized data
    w = taper(seq_out(n).T);
    seq_out(n).(datafield) ...
        = w' .* ((seq_in(n).(datafield) - mean_raw) ./ sqrt(var_raw));
end

% Get the current mean and variance of the tapered data
var_tapered = var([seq_out.(datafield)],0,2);  % Post-taper variance
mean_tapered = mean([seq_out.(datafield)],2);  % Post-taper mean
% Shift and rescale to conserve the raw mean and variance
for n = 1:length(seq_out)
    seq_out(n).(datafield) = sqrt(var_raw./var_tapered) .* (seq_out(n).(datafield) - mean_tapered) + mean_raw;
end
var_out = var([seq_out.(datafield)],0,2);  % Post-taper variance
mean_out = mean([seq_out.(datafield)],2);  % Post-taper mean
