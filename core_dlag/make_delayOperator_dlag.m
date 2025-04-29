function Q = make_delayOperator_dlag(params, f)
%
% Q = make_delayOperator_dlag(params, f)
%
% Constructs time-delay operator matrix for a given frequency.
%
% Arguments:
%
%     params -- DLAG model parameters
%     f      -- float; frequency, in units of (1/timeSteps)
%
% Outputs:
%
%     Q      -- GP time-delay operator matrix with dimensions 
%               (numGroups*xDim x xDim).
%
% Author: 
%     Evren Gokcen    egokcen@cmu.edu
% 
% Revision history:
%     07 Jun 2023 -- Initial full revision.

numGroups = length(params.xDim_within);
Q = cell(1,numGroups); % Time-delay operator matrix for each group

% Fill in Q_within (which is just the identity for each group)
Q_within = ones(sum(params.xDim_within),1);

% Fill in the rest of Q
for groupIdx = 1:numGroups
    % Fill in Q_across
    Q_across = [];
    if params.xDim_across > 0
       Q_across = exp(-1i*2*pi*f.*params.DelayMatrix(groupIdx,:)).';
    end
    % Concatenate across and within-group operators
    Q{groupIdx} = diag(vertcat(Q_across, Q_within));
end

% Collect time-delay operator matrices across groups
Q = vertcat(Q{:});
