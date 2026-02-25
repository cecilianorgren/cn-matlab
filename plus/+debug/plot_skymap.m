% Debug mms.plot_skymap
time = irf_time('2015-11-12T07:19:21.174999994Z','utc>epochtt')+1;
tind = find(abs(ePDist1.time-time)==min(abs(ePDist1.time-time)));
gseBref = mean(gseB1.tlim(time+[-0.005 0.005]).data);
energies = ePDist1(tind).depend{1};


nrows = 1;
ncols = 3;
npanels = ncols*nrows;

for isub = 1:npanels
  h(isub) = subplot(nrows,ncols,isub);
end

elevel = [12];
isub = 1;

plottype = 'sphere';
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,ePDist1.deflux,'tint',time,plottype,'energy',energies(elevel),'vectors',{gseBref,'B'});

plottype = 'flat';
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,ePDist1.deflux,'tint',time,plottype,'energy',energies(elevel),'vectors',{gseBref,'B'});

hca = h(isub); isub = isub + 1;
mms.plot_skymap_new(hca,ePDist1.deflux,'tint',time,plottype,'energy',energies(elevel),'vectors',{gseBref,'B'});  
