function blocks = mat2blocks(A, block_idxs)
% Helper function to convert block-diagonal matrix into cell array of
% blocks.
    numBlocks = length(block_idxs);
    blocks = cell(1, numBlocks);
    for blockIdx = 1:numBlocks
        currBlock = block_idxs{blockIdx};
        blocks{blockIdx} = A(currBlock(1):currBlock(2), currBlock(1):currBlock(2));
    end
end