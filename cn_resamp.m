% gsmE3
original=gsmE4;
nsamp=450;
freq=450*4;
rest=mod(size(original,1),freq);
n_ind=numel(1:size(original,1)-rest);
ind=1:size(original,1)-rest;

neworg_1=tocolumn(nanmean(reshape(original(ind,1),freq,n_ind/freq),1));
neworg_2=tocolumn(nanmean(reshape(original(ind,2),freq,n_ind/freq),1));
neworg_3=tocolumn(nanmean(reshape(original(ind,3),freq,n_ind/freq),1));
neworg_4=tocolumn(nanmean(reshape(original(ind,4),freq,n_ind/freq),1));

gsmE4lowres2=[neworg_1 neworg_2 neworg_3 neworg_4];