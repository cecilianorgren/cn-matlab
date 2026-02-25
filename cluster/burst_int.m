%%
data=load('BurstIntervals.txt');
%%
intervals=cell(size(data,1),2);
%% Storing burst intervals in cell array
for i=1:size(data,1) % rows
    for j=1:2 % columns
        intervals{i,j}=[data(i,1+6*(j-1)) data(i,2+6*(j-1))...
            data(i,3+6*(j-1)) data(i,4+6*(j-1)) data(i,5+6*(j-1)) data(i,6+6*(j-1))];
    end
end
%% Creating vector with 1 for burst mode and 0 for non burst mode
start_time=toepoch(intervals{1,1})-5;
stop_time=toepoch(intervals{end,end});
t_vec=start_time:stop_time;
vec=zeros(length(t_vec),2);
vec(:,1)=t_vec;
%%
for i=1:size(intervals,1)   
    t1=find(t_vec==toepoch(intervals{i,1}));
    t2=find(t_vec==toepoch(intervals{i,2}));
    vec(t1:t2,2)=1;
end

%%
bar(vec(1:10:end,1),vec(1:10:end,2),'y')
add_timeaxis;
%%
print -dpng burts_intervals.png