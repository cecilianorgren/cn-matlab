tint = irf.tint('2015-10-16T10:33:30.23Z',0.03);
sc = 1;
c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint).time).data);',sc); hatB0 = double(irf_norm(B0));
vectors = {hatB0,'B'};
%vectors(end+1,:) = {hatB0+[1 0 0],'B'};
%%
for ii=1:8;
  hca = subplot(2,4,ii);  
  mms.plot_skymap(hca,desDist1,'tint',tint,'energylevel',5+ii,'vectors',vectors);
  view(hca,[1 1 0])
end

%%
tint = irf.tint('2015-10-16T10:33:30.23Z',0.04);
mms.plot_projection(desDist1,'tint',tint,'xyz',[-1 0 0; 0 -1 0;hatB0]);

