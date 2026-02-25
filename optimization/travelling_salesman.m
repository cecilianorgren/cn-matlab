function [x_path,y_path] = travelling_salesman(stopsLon,stopsLat)
nStops = numel(stopsLon);
idxs = nchoosek(1:nStops,2);
%Calculate all the trip distances, assuming that the earth is flat in order to use the Pythagorean rule.

dist = hypot(stopsLat(idxs(:,1)) - stopsLat(idxs(:,2)), ...
             stopsLon(idxs(:,1)) - stopsLon(idxs(:,2)));
lendist = length(dist);

Aeq = spones(1:length(idxs)); % Adds up the number of trips
beq = nStops;

Aeq = [Aeq;spalloc(nStops,length(idxs),nStops*(nStops-1))]; % allocate a sparse matrix
for ii = 1:nStops
    whichIdxs = (idxs == ii); % find the trips that include stop ii
    whichIdxs = sparse(sum(whichIdxs,2)); % include trips where ii is at either end
    Aeq(ii+1,:) = whichIdxs'; % include in the constraint matrix
end
beq = [beq; 2*ones(nStops,1)];

intcon = 1:lendist;
lb = zeros(lendist,1);
ub = ones(lendist,1);

opts = optimoptions('intlinprog','Display','off');
[x_tsp,costopt,exitflag,output] = intlinprog(dist,intcon,[],[],Aeq,beq,lb,ub,opts);

segments = find(x_tsp); % Get indices of lines on optimal path

hca = gca;

x_path = stopsLon(idxs(find(x_tsp),2));
y_path = stopsLat(idxs(find(x_tsp),1));
hca = gca;
plot(hca,stopsLon(idxs(segments(1,:),:)),stopsLat(idxs(segments(1,:),:)),'-')
hold(hca,'on')
for iSeg = 2:numel(segments)
  plot(hca,stopsLon(idxs(segments(iSeg,:),:)),stopsLat(idxs(segments(iSeg,:),:)),'-')

end
hold(hca,'off')

[path_order] = construct_path(x_tsp,idxs);
return
%tours = detectSubtours(x_tsp,idxs);
numtours = length(tours); % number of subtours
fprintf('# of subtours: %d\n',numtours);
end
function [path_order] = construct_path(x_tsp,idxs)
  segments = find(x_tsp);
  nSegments = numel(segments)-1;
  path_order(1:2) = idxs(segments(1),:);
  indexSegments = idxs(segments(2:end),:);
  
  for iSeg = 1:nSegments    
    [I,J] = find(indexSegments==path_order(end));
    if isempty(I) % end of loop
      return
    end
    path_order(end+1) = indexSegments(I,setdiff([1 2],J));
    indexSegments = indexSegments(setdiff(1:size(indexSegments,1),I),:);
  end
  %path_order
    
end