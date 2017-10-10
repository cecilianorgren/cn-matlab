function out = squadd(in)
% Checks if timeseries has duplicate/overlapping time entries, and if so,
% adds the overlapping part to each other.

originalIn = in;
ttt = in(:,1);
%dt = diff(ttt);

out = [];
decreasing = 1;
while 1
    if ~isempty(in); dt = diff(in(:,1)); end
    decreasing = find(dt<0,1,'first');
    if isempty(decreasing); 
        out = [out; in]; % add last bit 
        break; 
    end
    ttt = in(:,1);
    te1 = decreasing;
    ts2 = decreasing+1;
    timee1 = ttt(te1); % end of first overlapping series
    times2 = ttt(ts2); % start of second overlapping series
    % find start of first, ts1
    %ts1 = find(abs(ttt(1:(te1-1))-times2)==min(abs(ttt(1:(te1-1))-times2)));
    %te2 = find(abs(ttt((ts2+1):end)-timee1)==min(abs(ttt((ts2+1):end)-timee1)));
    ts1 = find(ttt<ttt(ts2),1,'last');
    te2 = find(ttt>ttt(te1),1,'first');    
    [ts1 te1 ts2 te2]
    % added time series part    
    newNoOverlap = in(1:(ts1-1),:);
    newOverlap = irf_add(1,in(ts1:te1,:),1,in(ts2:te2,:));
    out = [out; newNoOverlap; newOverlap];  
    in = in(te2+1:end,:);
    %irf_plot({out,in},'comp')
    %irf_plot([out;in])
end
return;
for kk=1:numel(decreasing)
    te1 = decreasing(kk);
    ts2 = decreasing(kk)+1;
    timee1 = ttt(decreasing(kk)); % end of first overlappng series
    times2 = ttt(decreasing(kk)+1); % start of second overlapping series
    % find start of first, ts1
    ts1 = find(abs(ttt-times2)==min(abs(ttt-times2)));
    te2 = find(abs(ttt-timee1)==min(abs(ttt-timee1)));
    % added time series part
    new = irf_add(1,in(ts1:te1,:),1,in(ts2:te2,:));
    out(ts1:te2,:) = [];
    
end
    