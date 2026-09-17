function [P, uniqueGroups, iG_in_P] = allPartitions(N)

    % Generate all nonempty unique groups
    uniqueGroups = cell(2^N - 1, 1);

    for mask = 1:(2^N - 1)
        group = find(bitget(mask, 1:N));
        uniqueGroups{mask} = group;
    end

    % Generate all partitions recursively
    P = recursivePartitions(N);

    % Determine which unique-group indices belong to each partition
    iG_in_P = cell(size(P));

    for p = 1:numel(P)

        partition = P{p};
        groupIndices = zeros(1, numel(partition));

        for j = 1:numel(partition)

            group = sort(partition{j});

            % Convert group to its corresponding bitmask
            mask = sum(2.^(group - 1));

            % The mask is also the index in uniqueGroups
            groupIndices(j) = mask;
        end

        iG_in_P{p} = groupIndices;
    end
end


function P = recursivePartitions(N)

    if N == 1
        P = {{1}};
        return
    end

    Pprev = recursivePartitions(N-1);
    P = {};

    for i = 1:numel(Pprev)

        partition = Pprev{i};

        % Insert N into each existing group
        for j = 1:numel(partition)

            newPartition = partition;
            newPartition{j} = [newPartition{j}, N];

            P{end+1} = newPartition;
        end

        % Add N as a new singleton group
        P{end+1} = [partition, {[N]}];
    end
end

% function [P,groups,iG_in_P] = allPartitions(N)
% 
%   if N == 1
%     P = {{1}};
%     groups = {{1}};
%     return
%   end
% 
%   Pprev = allPartitions(N-1);
%   P = {};
%   iG_in_P = {};
% 
%   for i = 1:length(Pprev)
% 
%     partition = Pprev{i};
% 
%     % Put N into each existing group
%     for j = 1:length(partition)
%       newPartition = partition;
%       newPartition{j} = [newPartition{j}, N];
%       P{end+1} = newPartition;
% 
%       % Find unique-group index for every group
%       groupIndices = zeros(1, length(newPartition));
% 
%       for k = 1:length(newPartition)
%         groupIndices(k) = findGroupIndex( ...
%           newPartition{k}, uniqueGroups);
%       end
% 
%       iG_in_P{end+1} = groupIndices;
%     end
% 
%     % Or create a new group containing N
%     P{end+1} = [partition, {[N]}];
%   end
% 
%   groups = {};
% 
%   for i = 1:length(P)
%     for j = 1:length(P{i})
%       group = sort(P{i}{j});
%       key = mat2str(group);
% 
%       if ~any(strcmp(groups, key))
%           groups{end+1} = key;
%       end
%     end
%   end
%   groups = cellfun(@(x)str2num(x),groups,'UniformOutput',false);
% end
% % Helper function
% function iG = findGroupIndex(group, uniqueGroups)
% 
%     group = sort(group);
% 
%     iG = [];
% 
%     for i = 1:length(uniqueGroups)
% 
%         if isequal(group, sort(uniqueGroups{i}))
%             iG = i;
%             return
%         end
%     end
% 
%     error('Group not found in uniqueGroups.');
% end