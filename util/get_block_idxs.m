function block_idxs = get_block_idxs(group_sizes)
% Helper function to get indices that correpond to blocks of a matrix,
% according to group_sizes.
    numGroups = length(group_sizes);
    block_idxs = cell(1,numGroups);
    startIdx = 1; endIdx = 1;
    for groupIdx = 1:numGroups
        groupSize = group_sizes(groupIdx);
        endIdx = startIdx + groupSize - 1;
        block_idxs{groupIdx} = [startIdx endIdx];
        startIdx = endIdx + 1;
    end
end