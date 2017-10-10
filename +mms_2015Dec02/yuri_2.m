%% Set up figure
tintOverview = irf.tint('2015-12-02T01:14:53.00/2015-12-02T01:15:00.50Z');
tint = tintOverview;
%tintDist = irf.tint('2015-10-16T10:33:20.00Z',0.03);
[h1,h2] = initialize_combined_plot(10,3,2,3,'vertical');

%% Plot timeseries
if 0 % B, single sc
  tozoom = 1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?.tlim(tint).abs},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('B_{?,DMPA}',ic),'(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.99 0.9]);
else % B, 4 sc comp
  tozoom = 1:3;
  hca = irf_panel('Bx');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dmpaB1l2pre.x.tlim(tint),dmpaB2l2pre.x.tlim(tint),dmpaB3l2pre.x.tlim(tint),dmpaB4l2pre.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{x}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'sc1','sc2','sc3','sc4'},[0.99 0.9]);

  hca = irf_panel('By');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dmpaB1l2pre.y.tlim(tint),dmpaB2l2pre.y.tlim(tint),dmpaB3l2pre.y.tlim(tint),dmpaB4l2pre.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{y}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

  hca = irf_panel('Bz');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dmpaB1l2pre.z.tlim(tint),dmpaB2l2pre.z.tlim(tint),dmpaB3l2pre.z.tlim(tint),dmpaB4l2pre.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{z}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);
  
  hca = irf_panel('abs B');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dmpaB1l2pre.abs.tlim(tint),dmpaB2l2pre.abs.tlim(tint),dmpaB3l2pre.abs.tlim(tint),dmpaB4l2pre.abs.tlim(tint)},'comp');
  hca.YLabel.String = {'|B|','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);
end

  hca = irf_panel('vx');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ve1brst.x.tlim(tint),ve2brst.x.tlim(tint),ve3brst.x.tlim(tint),ve4brst.x.tlim(tint)},'comp');
  hca.YLabel.String = {'v_{x}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

  hca = irf_panel('vy');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ve1brst.y.tlim(tint),ve2brst.y.tlim(tint),ve3brst.y.tlim(tint),ve4brst.y.tlim(tint)},'comp');
  hca.YLabel.String = {'v_{y}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

  hca = irf_panel('vz');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ve1brst.z.tlim(tint),ve2brst.z.tlim(tint),ve3brst.z.tlim(tint),ve4brst.z.tlim(tint)},'comp');
  hca.YLabel.String = {'v_{z}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);
  
  hca = irf_panel('abs v');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ve1brst.abs.tlim(tint),ve2brst.abs.tlim(tint),ve3brst.abs.tlim(tint),ve4brst.abs.tlim(tint)},'comp');
  hca.YLabel.String = {'|v|','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

if 1
  hca = irf_panel('brst n');
  %c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1brst.tlim(tint),ne2brst.tlim(tint),ne3brst.tlim(tint),ne4brst.tlim(tint)},'comp');
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YScale = 'lin';
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);
end

if 1
  hca = irf_panel('sc Pot');
  %c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{P1brst.tlim(tint),P2brst.tlim(tint),P3brst.tlim(tint),P4brst.tlim(tint)},'comp');
  hca.YLabel.String = {'sc Pot','(V)'};
  hca.YScale = 'lin';
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);
end

if 0% Compare parallel and perpendicular phase space densities
clim = [-1 1];
ylim = [10 1e3];
for ic = 1:4
  if 1
    hca = irf_panel(irf_ssub('edistparperp ?',ic));
    c_eval('[~,hcb] = irf_spectrogram(hca,ePSDparperp?,''log'',''donotfitcolorbarlabel'');',ic)
    hold(hca,'on')  
    hold(hca,'off')
    irf_legend(hca,'(g)',[0.99 0.98],'color','k','fontsize',12)
    set(hca,'yscale','log');
    hca.CLim = clim;
    set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    hca.YLim = ylim;
    ylabel(hca,{'E_e','(eV)'});
    colormap(hca,cn.cmap('bluered3'))
    
    hhcb(isub) = hcb;
    if 1
      hcb.YTick = 0.6*[-1 1];
      hcb.YTickLabel = {'f_{perp}','f_{a/par}'};
      hcb.YLabel.String = irf_ssub('mms {?}',ic);
      hcb.YLabel.Color = mms_colors(irf_ssub('?',ic));
    end
    
    isub = isub + 1;    
  end
end
for ic = 1:4
  if 1
    hca = irf_panel(irf_ssub('edistparapar ?',ic));
    c_eval('[~,hcb] = irf_spectrogram(hca,ePSDparapar?,''log'',''donotfitcolorbarlabel'');',ic)
    irf_legend(hca,'(f)',[0.99 0.98],'color','k','fontsize',12)
    set(hca,'yscale','log');
    hca.CLim = clim;
    set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    hca.YLim = ylim;
    ylabel(hca,{'E','(eV)'});
    colormap(hca,cn.cmap('bluered3'))
    
    %hca.YLabel.String = {irf_ssub('E_{?}',ic),'(eV)'};
    
    % Colorbar labels
    hhcb(isub) = hcb;  
    if 0
      hcb.YLabel.Rotation = 0;
      hcb.YLabel.VerticalAlignment = 'middle';
      hcb.YLabel.String ={' ','f_{par}',' ','f_{apar}'};
      hcb.YLabel.Position = hcb.YLabel.Position + [1 0 0];
    end
    if 1
      hcb.YTick = 0.6*[-1 1];
      hcb.YTickLabel = {'f_{apar}','f_{par}'};
      hcb.YLabel.String = irf_ssub('mms {?}',ic);
      hcb.YLabel.Color = mms_colors(irf_ssub('?',ic));
    end
    
    isub = isub + 1;
  end
end
end

irf_zoom(h1(:),'x',tintOverview)
%irf_zoom(h1,'y')
irf_plot_axis_align

%% Plot single time particle distributions
ic = 2;

tintDist = irf.tint('2015-12-02T01:14:54.00Z/2015-12-02T01:14:59.00Z');
% tstep = 0.03*2;
% times = tintDist(1):tstep:tintDist(2);
c_eval('times = desDist?.tlim(tintDist).time;',ic)


for ii = 1:times.length;
  tint = times(ii)+[-1 1]*0.5*0.03*0.5;
  if exist('hmark'); delete(hmark); end
  hmark = irf_pl_mark(h1,tint.epochUnix','green');
  vlim = 15*1e3;
  elevlim = 20;
  strCMap = 'jet';
  %energies =  [30 220];
  projclim = [-1 4.5];
  palim = [1e-3 1e6];
  skymapEnergy = [65 278];
   
  c_eval('dist = desDist?;',ic)

  %c_eval('Vi0 = mean(vi?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  %hatVi0 = double(irf_norm(Vi0));
  c_eval('Ve0 = mean(ve?brst.resample(tint).data,1);',ic); 
  hatVe0 = double(irf_norm(Ve0));    
   
  % Get mean magnetic field direction
  c_eval('B0 = mean(gseB?l2pre.resample(dslE?brst.tlim(tint).time).data);',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);
  
  vectors = {hatB0,'B';hatE0,'E';hatVe0,'V_e'};%0;hatVe0,'V_e'};

  % Projection coordinate system
  x = hatB0;
  y = hatExB0;
  z = cross(x,y);
  
  isub = 1;
  
  % Plot psd 0 90 180
  hca = h2(isub); isub = isub + 1;
  psdtint = times(ii);tint;%+[-1 1]*0.03;
  c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?l2pre,''tint'',psdtint,''scPot'',P?brst,''ylim'',palim,''energies'',skymapEnergy);',ic)
  TitleStr = {irf_ssub('MMS ?',ic),[irf_time(psdtint.start.utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(psdtint.stop-psdtint.start) ' s']};
  %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];    
  %hca.Title.String = irf_ssub('C?',ic);
  hca.Title.String = TitleStr;
  
  % Plot skymap for a given energy
  hca = h2(isub); isub = isub + 1;      
  c_eval('mms.plot_skymap(hca,dist,''tint'',times(ii),''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
  %hca.Title.String = hca.Title.String{2};
  
  % Plot skymap for a given energy
  hca = h2(isub); isub = isub + 1;      
  c_eval('mms.plot_skymap(hca,dist,''tint'',times(ii),''energy'',skymapEnergy(2),''vectors'',vectors,''flat'');',ic)
  %hca.Title.String = hca.Title.String{2};

  % Plot project ion onto a plane
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',times(ii),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)

  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',times(ii),'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)
  
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',times(ii),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)
    
  pause(0.1)
  cn.print([irf_ssub('BvnP_psds_mms?_',ic) irf_time(times(ii),'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
end