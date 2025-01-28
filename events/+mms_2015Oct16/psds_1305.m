%% Load data

ic = 1:4;
tintOverview = irf.tint('2015-10-16T13:05:24.00Z/2015-10-16T13:05:47.00Z'); % magnetosphere-magnetosheath-magnetosphere
tintPsd = irf.tint('2015-10-16T13:05:38.00Z',0.03); tint = tint(1);
tint = tintOverview;
if 0
c_eval('Psp?brst=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_psp'',tint);',ic); % prope to spacecraft potential
c_eval('tic; desDist? = mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint); toc',ic);
c_eval('dmpaB?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',ic);
c_eval('dmpaB?brst=mms.db_get_ts(''mms?_dfg_brst_ql'',''mms?_dfg_brst_dmpa'',tint);',ic)
c_eval('vex?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkX'',tint);',ic);
c_eval('vey?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkY'',tint);',ic);
c_eval('vez?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkZ'',tint);',ic);
c_eval('ve?brst=irf.ts_vec_xyz(vex?brst.time,[vex?brst.data vey?brst.data vez?brst.data]);',ic)
c_eval('ve?brst.name = ''mms? ve brst'';')
end

%% Initialize plot
if 0
  nrows = 3;
  ncols = 3;
  h1 = subplot(nrows,ncols,[1 2 3]);

  for ii = 4:9;
    h2(ii) = subplot(nrows,ncols,ii);
  end
else
  [h1,h2] = mms.initialize_comined_plot(5,3,2,3,'vertical');
end

%% Plot time series
ic = 3;
tint = tintOverview;
if 1 % B, single sc
  hca = irf_panel(irf_ssub('B?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);
end

if 1 % Omnidirectional   
  hca = irf_panel(irf_ssub('DEF omni ?',ic));
  c_eval('[~,eomnicb] = irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hold(hca,'on');
  c_eval('irf_plot(hca,Psp?brst.tlim(tint)*(-1),''color'',''black'');',ic)  
  hold(hca,'off');
  %irf_legend(hca,'(e)',[0.99 0.98],'color','black','fontsize',12)
  irf_legend(hca,{'Psp'},[0.01 0.4],'color','black')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %caxis(hca,[-15 -5])
  hca.YLim = [10 1e3];
  ylabel(hca,{'E_{e,OMNI}',irf_ssub('(eV) mms ?',ic)},'Interpreter','tex');
  %eomnicb.YLabel.String = 'DEF';
  colormap(hca,'jet')
  %hca.CLim = [4.5 9];
end

if 1 % parallel DEF   
  hca = irf_panel(irf_ssub('DEF par ?',ic));
  c_eval('[~,eparcb] = irf_spectrogram(hca,eDEFpar?,''log'',''donotfitcolorbarlabel'');',ic)
  hold(hca,'on');
  c_eval('irf_plot(hca,Psp?brst.tlim(tint)*(-1),''color'',''black'');',ic)  
  hold(hca,'off');
  %irf_legend(hca,'(e)',[0.99 0.98],'color','black','fontsize',12)
  irf_legend(hca,{'Psp'},[0.01 0.4],'color','black')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %caxis(hca,[-15 -5])
  hca.YLim = [10 1e3];
  ylabel(hca,{'E_{e,par}',irf_ssub('(eV) mms ?',ic)},'Interpreter','tex');
  %eomnicb.YLabel.String = 'DEF';
  colormap(hca,'jet')
  %hca.CLim = [4.5 9];
end

if 1 % antiparallel DEF
  hca = irf_panel(irf_ssub('DEF apar ?',ic));
  c_eval('[~,eomnicb] = irf_spectrogram(hca,eDEFapar?,''log'',''donotfitcolorbarlabel'');',ic)
  hold(hca,'on');
  c_eval('irf_plot(hca,Psp?brst.tlim(tint)*(-1),''color'',''black'');',ic)  
  hold(hca,'off');
  %irf_legend(hca,'(e)',[0.99 0.98],'color','black','fontsize',12)
  irf_legend(hca,{'Psp'},[0.01 0.4],'color','black')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %caxis(hca,[-15 -5])
  hca.YLim = [10 1e3];
  ylabel(hca,{'E_{e,apar}',irf_ssub('(eV) mms ?',ic)},'Interpreter','tex');
  %eomnicb.YLabel.String = 'DEF';
  colormap(hca,'jet')
  %hca.CLim = [4.5 9];
end

if 1 % perpendicular DEF
  hca = irf_panel(irf_ssub('DEF perp ?',ic));
  c_eval('[~,eomnicb] = irf_spectrogram(hca,eDEFperp?,''log'',''donotfitcolorbarlabel'');',ic)
  hold(hca,'on');
  c_eval('irf_plot(hca,Psp?brst.tlim(tint)*(-1),''color'',''black'');',ic)  
  hold(hca,'off');
  %irf_legend(hca,'(e)',[0.99 0.98],'color','black','fontsize',12)
  irf_legend(hca,{'Psp'},[0.01 0.4],'color','black')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  %caxis(hca,[-15 -5])
  hca.YLim = [10 1e3];
  ylabel(hca,{'E_{e,perp}',irf_ssub('(eV) mms ?',ic)},'Interpreter','tex');
  %eomnicb.YLabel.String = 'DEF';
  colormap(hca,'jet')
  %hca.CLim = [4.5 9];
end


irf_zoom(h1,'x',tintOverview)
irf_zoom(h1(1),'y')
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);


%% Plot particle distributions

tint = tintPsd+0.1*[-1 1]+0*(-1.2)+20/8*0.3;

if exist('hmark'); delete(hmark); end
hmark = irf_pl_mark(h1,tint.epochUnix','black');
vlim = 10*1e3; 
elevlim = 20;
strCMap = 'jet';
%energies =  [30 220];
projclim = [-2 5.2];
palim = [1e-3 1e6];


c_eval('dist = desDist?;',ic)

%c_eval('Vi0 = mean(vi?brst.resample(dmpaB?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
%hatVi0 = double(irf_no rm(Vi0));
%c_eval('Ve0 = mean(ve?brst.resample(dmpaB?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
%hatVe0 = double(irf_norm(Ve0));

c_eval('PsP = mean(Psp?brst.tlim(tint).data);',ic); 
skymapEnergy = [-PsP*0.75 -PsP*1.1];

% Get mean magnetic field direction
c_eval('B0 = mean(dmpaB?brst.tlim(tint+0.05*[-1 1]).data);',ic); 
hatB0 = double(irf_norm(B0));   
vectors = {hatB0,'B'};%[1 0 0],'X'};%0;hatVe0,'V_e'};



z = hatB0;
x = cross(cross(z,[1 0 0]),z); x = x/norm(x);  
y = cross(z,x); y = y/norm(y);
x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

isub = 1;

% Plot psd 0 90 180
hca = h2(isub); isub = isub + 1;
c_eval('mms.plot_pitchangles(hca,dist,dmpaB?brst,''tint'',tint,''scPot'',Psp?brst,''ylim'',palim,''energies'',skymapEnergy);',ic)
hca.Title.String = {irf_ssub('MMS ?',ic), [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s']};    

% Plot skymap for a given energy
hca = h2(isub); isub = isub + 1;      
c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
hca.Title.String = hca.Title.String{2};
%hca.Title.String = '';

% Plot skymap for a given energy
hca = h2(isub); isub = isub + 1;      
c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(2),''vectors'',vectors,''flat'');',ic)
hca.Title.String = hca.Title.String{2};
%hca.Title.String = '';

% Plot project ion onto a plane
hca = h2(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim,'vlabel',{'x','y','z'});
hca.Title.String = '';
colormap(hca,strCMap)

hca = h2(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim,'vlabel',{'z','x','y'});
hca.Title.String = '';
colormap(hca,strCMap)

hca = h2(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim,'vlabel',{'y','z','x'});
hca.Title.String = '';
colormap(hca,strCMap)

% Small adjustments
for ii= 2:5; 
  grid(h1(ii),'off'); 
  h1(ii).CLim = [4.5 8.5];  
end
h2(2).CLim = h2(2).CLim*0.1;
for ii = 4:6
  h2(ii).CLim = [0 5];
end

%h2(3).CLim = [0 2e3];
%pause(1)
%cn.print([irf_ssub('mms?_',ic) irf_time(tint(1),'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
