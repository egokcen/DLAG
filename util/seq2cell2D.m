function dat = seq2cell2D(seq, dims, varargin)
%
% dat = seq2cell2D(seq, dims, ...) 
%
% Description: Convert data in sequence format to a cell array of data 
%              matrices. Each cell array entry corresponds to a group,
%              and data is concatenated across trials.
%
% Arguments:
%
%     Required:
%
%     seq    -- structure whose nth entry (corresponding to the nth sequence)
%               has (relevant) fields
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
%     dat    -- (1 x numGroups) cell array; list of data matrices 
%              {(dims(1) x numSamples), (dims(2) x numSamples), ...}
%
% Authors: 
%     Evren Gokcen    egokcen@cmu.edu

datafield = 'y';
extraOpts = assignopts(who, varargin);
numGroups = length(dims);

group_blocks = get_block_idxs(dims);
dat = cell(1,numGroups);
N = length(seq);
for groupIdx = 1:numGroups
    currGroup = group_blocks{groupIdx}(1):group_blocks{groupIdx}(2);
    dat{groupIdx} = [];
    for n = 1:N
        dat{groupIdx} = [dat{groupIdx} seq(n).(datafield)(currGroup,:)];
    end
end