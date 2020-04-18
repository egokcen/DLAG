function timeIdxs = time_shift(T, delay, max_delay)
% Description: Given a time window T (in ms), select an interval between 1
% and T, based on delay, of length T - 2*max_delay.
   
    % Base time window (i.e., 0 time shift)
    startIdx = max_delay + 1;
    endIdx = T - max_delay;
    % Negative delay means shift forward in time
    startIdx = startIdx - delay;
    endIdx = endIdx - delay;
    timeIdxs = startIdx:endIdx;
end