function numWorkers = StartParPool(numParRuns)

% Helper function to configure MATLAB's parallelization tools
% (e.g., parfor)

    currParPool = gcp('NoCreate');

    numWorkers = min(feature('NumCores'), numParRuns);

    if isempty(currParPool)
        parpool( numWorkers );
    elseif currParPool.NumWorkers < numWorkers
        delete(currParPool)
        parpool( numWorkers );
    end

end
