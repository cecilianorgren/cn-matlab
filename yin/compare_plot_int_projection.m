%% Plot overview figure with focus on electrons, including single time electron distributions
npanels = 3;
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tintZoom = tint;

cmap = 'jet';
[h1,h2] = initialize_combined_plot(npanels,3,3,0.4,'vertical');
ic = 1;
iisub = 0;
cmap = colormap('jet');

isub = 0;
zoomy = [];

if 1 % B  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','x','|B|'},[0.98 0.9],'fontsize',12);  
end
if 1 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 0 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryE?.x,bdryE?.y,bdryE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);  
end
if 0 % E perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x,gseE?perp.y,gseE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % J  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
end
if 0 % Vi  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVi?.x.tlim(tint),bdryVi?.y.tlim(tint),bdryVi?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?.x.tlim(tint),bdryVe?.y.tlim(tint),bdryVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);     
end
if 0 % ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);   
end
if 0 % Ve par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(km/s)'};
  %hca.YLim = [-1100 1100];  
end
if 0 % Ve perp par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?perp.x.tlim(tint),bdryVe?perp.y.tlim(tint),bdryVe?perp.z.tlim(tint),bdryVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp}','n_{\perp}','B_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % e DEF omni 64
  isub = isub + 1;
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
if 0 % ePDist pa 64
  isub = isub + 1;
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [100 30000];  
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % Te par perp Ti/Tref
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h1(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h1(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h1,'x',tintZoom)
irf_zoom(h1(zoomy),'y')
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);
if 0
hmark = irf_pl_mark(h1(1:6),tintBCS, 'yellow');
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
end

%% Plot single time particle distributions, 1 sc, 4 projections,
tintDist = irf.tint('2017-07-06T13:54:00.00Z/2017-07-06T13:54:40.00Z'); 
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('dist_scm = ePDist?;',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
%c_eval('ePitch = ePitch?.convertto(''s^3/km^6'');',ic)

% Plot format input
vlim = 50*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [-3 1.5];  
projclim_int = [-13 -9];
  

c_eval('times = ePDist?.time;',ic)
[tind,~] = times.tlim(tintDist);
clear hmark_tmp hmark
doReducedF = 1;
for it_ = 1:30:numel(tind);
  it = tind(it_);
  time = times(it);
  it
  if exist('hmark'); delete(hmark); end
  c_eval('hmark_tmp = irf_pl_mark(h1(?),time.epochUnix,''green''); hmark(?) = hmark_tmp;',1:npanels)
  
%   c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
%   c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
%   c_eval('hatExB = cross(hatE,hatB);',ic)
%   par = hatB;
%   perp1 = hatExB;
%   perp2 = cross(par,perp1);  
  
  c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
  par = hatB;
  newy = [0 1 0]; 
  perp1 = cross(newy,par);
  perp2 = cross(par,perp1);  
  

  timeUTC = time.utc;      
  isub = 1;

  %vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B'};
  vectors = {[1 0 0],'x'; [0 1 0],'y'; [0 0 1],'z'};
  
  if 1 % Perpendicular plane, slice
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{Bx[0 1 0]}','v_{Bx([0 1 0]xB)}','v_{B}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'vlabel',vlabels);
    hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
  end 
  if 1 % B plane 1, slice
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{Bx[0 1 0]}','v_{B}','-v_{Bx([0 1 0]xB)}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'vlabel',vlabels);        
    hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
  end 
  if 1 % B plane 2, slice
    hca = h2(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{Bx([0 1 0]xB)}','v_{B}','v_{Bx[0 1 0]}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'vlabel',vlabels);        
    hca.Title.String = sprintf('theta_{elev} = %g deg',elevlim);
  end
  vint = [-Inf Inf];
  if 1 % Perpendicular plane, integrated
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{Bx[0 1 0]}','v_{Bx([0 1 0]xB)}','v_{B}'};
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end 
  if 1 % B plane 1, integrated
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{Bx[0 1 0]}','v_{B}','-v_{Bx([0 1 0]xB)}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end  
  if 1 % B plane 2, integrated
    hca = h2(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{Bx([0 1 0]xB)}','v_{B}','v_{Bx[0 1 0]}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end  
  vint = 10000*[-1 1];
  if 1 % Perpendicular plane, integrated, smaller vint
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{Bx[0 1 0]}','v_{Bx([0 1 0]xB)}','v_{B}'};
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end 
  if 1 % B plane 1, integrated, smaller vint
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{Bx[0 1 0]}','v_{B}','-v_{Bx([0 1 0]xB)}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end  
  if 1 % B plane 2, integrated, smaller vint
    hca = h2(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{Bx([0 1 0]xB)}','v_{B}','v_{Bx[0 1 0]}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = sprintf('vint = [%g %g] km/s',vint(1),vint(2));
  end  
  
  if 0 % Pitchangle distribution
    hca = h2(isub); isub = isub + 1;
    plot(hca,ePitch.depend{1}(it,:),squeeze(ePitch.data(it,:,[1 ceil(numel(ePitch.depend{2})/2) numel(ePitch.depend{2})])));
    hca.YScale = 'log'; hca.XScale = 'log';
    hca.YLabel.String = ['f_e (' ePitch.units ')'];
    hca.XLabel.String = 'E (eV)';
    hca.XLim = [50 30000];
    legend(hca,{'0','90','180'})
    hca.YTick = 10.^[-4:2];
    hca.YLim = [1e-4 1e2];
  end
  if 0 % invisible
    hca = h2(isub); isub = isub + 1;
    hca.Visible = 'off';
  end
  if 0 % invisible
    hca = h2(isub); isub = isub + 1;
    hca.Visible = 'off';
  end
  % ExB plane, with markings
  if 0
    vectors = {hatExB,'ExB'; hatE,'E';hatB,'B';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z'};
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{E,\perp}','v_B'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    hca.Title.String = '';
  end
  
  h2(1).Title.String = {timeUTC(1:23),h2(1).Title.String};
  for ii = 1:4
    colormap(h2(ii),strCMap)
  end
  %cn.print(['e_proj_in_fix_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','path',[eventPath 'proj_int_mms1/'])
end
