function block_mask = create_block_mask(group_sizes)
% Helper function to create mask to keep only blocks of a matrix.
% Used in projected gradient step.
    numGroups = length(group_sizes);
    blocks = cell(1,numGroups);
    for groupIdx = 1:numGroups
        blocks{groupIdx} = ones(group_sizes(groupIdx)); 
    end
    block_mask = blkdiag(blocks{:});
end