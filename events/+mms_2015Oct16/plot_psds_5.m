tintOverview = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tint = tintOverview;
%tintDist = irf.tint('2015-10-16T10:33:20.00Z',0.03);
[h1,h2] = mms.initialize_comined_plot(9,3,2,3,'vertical');

ic  = 1;

if 1 % B, single sc
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);
end
if 1
  hca = irf_panel('brst E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'E_x','E_y','E_z'},[0.95 0.95]);
end
if 0
  %%
  Bfilt=irf_filt(Bscm,100,4000,8192,3);
  hca = irf_panel('brst B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Bfilt.tlim(tint).x,Bfilt.tlim(tint).y},''comp'');',ic)
  hca.YLabel.String = {'B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z'},[0.95 0.95]);
  
end
if 1  
  hca = irf_panel('brst ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{ve?brst.tlim(tint).x,ve?brst.tlim(tint).y,ve?brst.tlim(tint).z},''comp'');',ic);
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);
end

if 1  
  hca = irf_panel('brst vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{vi?_lowres.tlim(tint).x,vi?_lowres.tlim(tint).y,vi?_lowres.tlim(tint).z},''comp'');',ic);
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);
end
if 1
  hca = irf_panel('brst n');
  %c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,ne?_lowres.tlim(tint));',ic);
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YScale = 'lin';
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'n_e'},[0.95 0.95]);
end

if 1
  hca = irf_panel('brst Te');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Te?_lowres.tlim(tint).x,Te?_lowres.tlim(tint).y,Te?_lowres.tlim(tint).z},''comp'');',ic);
  hca.YLabel.String = {'T','(eV)'};
  hca.YScale = 'lin';
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'T_{x}','T_{y}','T_{z}'},[0.95 0.95]);
end

%c
clim = [-1 1];
ylim = [10 1e3];

if 1
  
    hca = irf_panel(irf_ssub('edistparperp ?',ic));
    c_eval('[~,hcb] = irf_spectrogram(hca,ePSDparperp?,''log'',''donotfitcolorbarlabel'');',ic)
    hold(hca,'on')  
    hold(hca,'off')
    %irf_legend(hca,'(g)',[0.99 0.98],'color','k','fontsize',12)
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
end
if 1
  hca = irf_panel(irf_ssub('edistparapar ?',ic));
  c_eval('[~,hcb] = irf_spectrogram(hca,ePSDparapar?,''log'',''donotfitcolorbarlabel'');',ic)
  %irf_legend(hca,'(f)',[0.99 0.98],'color','k','fontsize',12)
  set(hca,'yscale','log');
  hca.CLim = clim;
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLim = ylim;
  ylabel(hca,{'E','(eV)'});
  colormap(hca,cn.cmap('bluered3'))

  %hca.YLabel.String = {irf_ssub('E_{?}',ic),'(eV)'};

  % Colorbar labels
  if 1
    hcb.YTick = 0.6*[-1 1];
    hcb.YTickLabel = {'f_{apar}','f_{par}'};
    hcb.YLabel.String = irf_ssub('mms {?}',ic);
    hcb.YLabel.Color = mms_colors(irf_ssub('?',ic));
  end
end
  
if 0
  hca = irf_panel('brst Te par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Te?par},''comp'');',ic);
  hca.YLabel.String = {'T','(eV)'};
  hca.YScale = 'lin';
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'T_{||}','T_{y}','T_{z}'},[0.95 0.95]);
end



irf_zoom(h1,'x',tintOverview)
irf_zoom(h1(1:5),'y')
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic)

%% Plot single time particle distributions
tintDist = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:35.00Z');
tints = tintDist(1) + [1 2 3 4 5 6 7:0.03:12 13 14 15 16:0.2:21 22:0.2:23.8 24:0.03:29.3 30:1:38];
tint = tintDist(1)+0.1;
ii=0;
%while tint<tintDist(end)
for ii = 1;117;40;17:48;tints.length;
  %%
  tint = irf.tint('2015-10-16T10:33:29.40Z',1);
  tint = tint(1)
  %ii= ii+1;
  %if ii>10; break; end
  %tint = tint + 0.1; 
  %tint = tint(ii);
  if exist('hmark'); delete(hmark); end
  hmark = irf_pl_mark(h1,tint.epochUnix','green');
  vlim = 20*1e3;
  elevlim = 20;
  strCMap = 'jet';
  %energies =  [30 220];
  projclim = [-2 5.2];
  palim = [1e-3 1e6];
  skymapEnergy = [50 138];
  
  c_eval('dist = desDist?;',ic)

  c_eval('Vi0 = mean(vi?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  hatVi0 = double(irf_norm(Vi0));
  c_eval('Ve0 = mean(ve?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  hatVe0 = double(irf_norm(Ve0));

  % Get mean magnetic field direction
  c_eval('B0 = mean(dmpaB?brst.resample(dslE?brst.tlim(tint+[-0.01 0.01]).time).data);',ic); 
  hatB0 = double(irf_norm(B0));   
  vectors = {hatB0,'B';[0 0 1],'Z'};%0;hatVe0,'V_e'};
  
   
  z = hatB0;
  x = cross(cross(z,[1 0 0]),z); x = x/norm(x);  
  y = cross(z,x); y = y/norm(y);
  %mva
  
  isub = 1;
  
  % Plot psd 0 90 180
  hca = h2(isub); isub = isub + 1;
  c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?brst,''tint'',tint+0.03*[-1 1],''scPot'',P?,''ylim'',palim,''energies'',skymapEnergy);',ic)
  %hca.Title.String = [irf_time(tint(1).utc,'utc>utc_yyyy-mm-ddTHH:MM:SS.mmm') ' + ' num2str(tint.stop-tint.start) ' s'];    

  % Plot skymap for a given energy
  hca = h2(isub); isub = isub + 1;      
  c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(1),''vectors'',vectors,''flat'');',ic)
  %hca.Title.String = hca.Title.String{2};
  
  % Plot skymap for a given energy
  hca = h2(isub); isub = isub + 1;      
  c_eval('mms.plot_skymap(hca,dist,''tint'',tint,''energy'',skymapEnergy(2),''vectors'',vectors,''flat'');',ic)
  %hca.Title.String = hca.Title.String{2};

  % Plot project ion onto a plane
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)

  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)
  
  hca = h2(isub); isub = isub + 1; 
  mms.plot_projection(hca,dist,'tint',tint(1),'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'vectors',vectors,'clim',projclim);
  hca.Title.String = '';
  colormap(hca,strCMap)
  
  grid(h1(6),'off')
  pause(1)
  %cn.print([irf_ssub('mms?_',ic) irf_time(tint(1),'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
end

if 0
for ii = 1:numel(h2)
  colormap(h2(ii),'jet')
  h2(ii).CLim = [-2 5];
end
end 