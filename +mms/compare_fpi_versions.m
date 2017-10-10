cmap = irf_colormap('poynting');
%%
tint = irf.tint('2015-10-16T10:33:47.40Z',0.03); % magnetosphere-magnetosheath-magnetosphere


nTimes = 4;
energy = 90;
ncols = nTimes;
nrows = 3;
dt = 0.03;
%clim = [0 23000];
time = tint(1);

isub = 1;
for it = 1:nTimes  
  
  time = tint(1) + it*dt;
  hca = subplot(nrows,ncols,isub); h(isub)= hca; isub = isub + 1;  
  mms.plot_skymap(hca,desDist,'tint',time,'energy',energy,'flat'); 
  %view([0 0 1])    
  %hca.CLim = clim;
  hca.Title.String = {desDist.userData.GlobalAttributes.Logical_file_id{1}(end-5:end), hca.Title.String{:}};
end
for it = 1:nTimes 
  time = tint(1) + it*dt;
  hca = subplot(nrows,ncols,isub); h(isub)= hca; isub = isub + 1;  
  mms.plot_skymap(hca,desDist1,'tint',time,'energy',energy,'flat'); 
  %view([0 0 1])  
  %hca.CLim = clim;
  hca.Title.String = {desDist1.userData.GlobalAttributes.Logical_file_id{1}(end-5:end), hca.Title.String{:}};
end
for it = 1:nTimes 
  time = tint(1) + it*dt;
  hca = subplot(nrows,ncols,isub); isub = isub + 1;  
  mms.plot_skymap(hca,(desDist-desDist1),'tint',time,'energy',energy,'flat'); 
  %view([0 0 1])  
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colormap(hca,cmap)
  hca.Title.String = {[desDist.userData.GlobalAttributes.Logical_file_id{1}(end-5:end) '-' desDist1.userData.GlobalAttributes.Logical_file_id{1}(end-5:end)], hca.Title.String{:}};
end
cmax = 0;
for ii = 1:8; 
  tmpcmax = max(h(ii).CLim);
  if tmpcmax>cmax
    cmax = tmpcmax;
  end
end
for ii = 1:8; h(ii).CLim = [0 cmax]; end

%cn.print(['des-dist_' [desDist.userData.GlobalAttributes.Logical_file_id{1}(end-5:end) '_' desDist1.userData.GlobalAttributes.Logical_file_id{1}(end-5:end)] '_' irf_time(tint(1),'epochtt>utc_yyyymmddTHHMMSS') '_approx' num2str(energy) 'eV']);
%%
tint = irf.tint('2015-10-16T10:33:28.000Z',0.03); % magnetosphere-magnetosheath-magnetosphere
time = tint(1);
isub = 1;
for it = 1:(nrows*ncols)
  time = tint(1) + it*dt;
  hca = subplot(nrows,ncols,isub); isub = isub + 1;  
  mms.plot_skymap(hca,(desDist-desDist1),'tint',time,'energy',energy,'flat'); 
  %view([0 0 1])  
  
  hca.Title.String = {[desDist.userData.GlobalAttributes.Logical_file_id{1}(end-5:end) '-' desDist1.userData.GlobalAttributes.Logical_file_id{1}(end-5:end)], hca.Title.String{:}};
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
irf_colormap('poynting')


%%
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:38.00Z');
ic = 4;

h = irf_plot(7);
isub = 1;

hca = h(isub); isub = isub + 1; 
c_eval('irf_plot(hca,{ve?brst.x,ve?brst_new.x},''comp'');',ic)
irf_legend(hca,{'old','new'},[0.02 0.98])
hca.YLabel.String = {'v_{e,x}','[km/s]'};

hca = h(isub); isub = isub + 1; 
c_eval('irf_plot(hca,{ve?brst.y,ve?brst_new.y},''comp'');',ic)
irf_legend(hca,{'old','new'},[0.02 0.98])
hca.YLabel.String = {'v_{e,y}','[km/s]'};

hca = h(isub); isub = isub + 1; 
c_eval('irf_plot(hca,{ve?brst.z,ve?brst_new.z},''comp'');',ic)
irf_legend(hca,{'old','new'},[0.02 0.98])
hca.YLabel.String = {'v_{e,z}','[km/s]'};

hca = h(isub); isub = isub + 1; 
c_eval('irf_plot(hca,{ne?brst,ne?brst_new},''comp'');',ic)
irf_legend(hca,{'old','new'},[0.02 0.98])
hca.YLabel.String = {'n_{e}','[cm^{-3}]'};

if 1 
  hca = h(isub); isub = isub + 1; 
  c_eval('irf_plot(hca,{Te?brst.x,Te?brst_new.x},''comp'');',ic)
  irf_legend(hca,{'old','new'},[0.02 0.98])
  hca.YLabel.String = {'T_{e,xx}','[eV]'};

  hca = h(isub); isub = isub + 1; 
  c_eval('irf_plot(hca,{Te?brst.y,Te?brst_new.y},''comp'');',ic)
  irf_legend(hca,{'old','new'},[0.02 0.98])
  hca.YLabel.String = {'T_{e,yy}','[eV]'};

  hca = h(isub); isub = isub + 1; 
  c_eval('irf_plot(hca,{Te?brst.z,Te?brst_new.z},''comp'');',ic)
  irf_legend(hca,{'old','new'},[0.02 0.98])
  hca.YLabel.String = {'T_{e,zz}','[eV]'};
end
if 0
  hca = h(isub); isub = isub + 1; 
  c_eval('irf_plot(hca,{Pe?brst.x*1e8,Pe?brst_new.x},''comp'');',ic)
  irf_legend(hca,{'old','new'},[0.02 0.98])
  hca.YLabel.String = {'P_{e,xx}','[erg/cm^3]'};

  hca = h(isub); isub = isub + 1; 
  c_eval('irf_plot(hca,{Pe?brst.y*1e8,Pe?brst_new.y},''comp'');',ic)
  irf_legend(hca,{'old','new'},[0.02 0.98])
  hca.YLabel.String = {'P_{e,yy}','[erg/cm^3]'};

  hca = h(isub); isub = isub + 1; 
  c_eval('irf_plot(hca,{Pe?brst.z*1e8,Pe?brst_new.z},''comp'');',ic)
  irf_legend(hca,{'old','new'},[0.02 0.98])
  hca.YLabel.String = {'P_{e,zz}','[erg/cm^3]'};
end

h(1).Title.String = irf_ssub('MMS ?',ic);

irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')

%%
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:38.00Z');
ic = 4;

h = irf_plot(3);
isub = 1;

hca = h(isub); isub = isub + 1; 
irf_plot(hca,'ne?brst_new','comp');
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.02 0.98])
hca.YLabel.String = {'n_{e,new}','[cm^{-3}]'};

hca = h(isub); isub = isub + 1; 
irf_plot(hca,'ne?brst','comp');
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.02 0.98])
hca.YLabel.String = {'n_{e,old}','[cm^{-3}]'}; 

hca = h(isub); isub = isub + 1; 
irf_plot(hca,'ne?_lowres','comp');
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.02 0.98])
hca.YLabel.String = {'n_{e,old,filtered}','[cm^{-3}]'}


%%
% Gradient of electron pressure
% start of this time interval is a relatively quiet period
tref = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:40.00Z')+[1 0]; 
c_eval('Pref? = Pe?_lowres.tlim(tref).data(1,:);') % Pa
c_eval('Pe?rel1 = Pe?_lowres+Pref1-Pref?;') % Remove offset
c_eval('r? = [Pe1rel1.time.epochUnix gseR?.resample(Pe1rel1.time).data Pe?rel1.abs.resample(Pe1rel1.time).data/3];')
for ii = 1:size(r1,1)
  gradPhie(ii,1)=r1(ii,1);
  gradPhie(ii,2:4)=c_4_gradphi(r1(ii,:),r2(ii,:),r3(ii,:),r4(ii,:));
end
gradPe = irf.ts_vec_xyz(irf_time(gradPhie(:,1),'epoch>utc'),gradPhie(:,2:4));
gradPe.units = 'Pa/km';
gradPe.name = 'grad(Pe)';