function [lengths, idxCells] = clusterLengths(x)

    % Ensure row vector
    x = x(:)';

    % Detect transitions
    dx = diff([0 x 0]);

    % Start and end of clusters of 1s
    startIdx = find(dx == 1);
    endIdx   = find(dx == -1) - 1;

    % Lengths of clusters
    lengths = endIdx - startIdx + 1;

    % Build cell array with indices of each cluster
    nClusters = numel(startIdx);
    idxCells = cell(1,nClusters);

    for k = 1:nClusters
        idxCells{k} = startIdx(k):endIdx(k);
    end

end