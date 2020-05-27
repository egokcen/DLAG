function groupSeq = partitionObs(seq, dims, varargin)
%
% groupSeq = partitionObs(seq, dims, varargin)
%
% Description: Partition sequences of observations according to the groups 
%              given in dims.
%
% Arguments:
%
%     Required:
%
%     seq    -- data structure, whose nth entry (corresponding to
%               the nth trial) has fields
%                   trialId               -- unique trial identifier
%                   T (1 x 1)             -- number of timesteps
%                   (datafield) (dim x T) -- continuous valued data
%
%     dims   -- (1 x numGroups) array; dimensionalities of each group
%
%     Optional:
%     
%     datafield -- string; Name of data field in seq (default: 'y')
%
% Outputs:
%
% Authors:
%     Evren Gokcen    egokcen@cmu.edu
%
% Revision history:
%     16 May 2020 -- Initial full revision.

datafield = 'y';
extraOpts = assignopts(who, varargin);

numGroups = length(dims);
group_blocks = get_block_idxs(dims);
groupSeq = cell(1,numGroups);
N = length(seq);
for groupIdx = 1:numGroups
    currGroup = group_blocks{groupIdx}(1):group_blocks{groupIdx}(2);
    for n = 1:N
        groupSeq{groupIdx}(n).trialId = seq(n).trialId;
        groupSeq{groupIdx}(n).T = seq(n).T;
        groupSeq{groupIdx}(n).(datafield) = seq(n).(datafield)(currGroup,:); 
    end
end