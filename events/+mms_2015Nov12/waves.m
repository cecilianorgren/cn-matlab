%% Waves at parallel electric field in current sheet
%% Minimum variance anlaysis
% GSE
v1 = [0.7460 0.3689 0.5545];
v2 = [0.6552 -0.5560 -0.5115];
v3 = [0.1196 0.7449 -0.6564];

%% Spectrogram
c_eval('wavE?hmfe = irf_wavelet(dslE?hmfe.abs.tlim(tintZoom),''wavelet_width'',5.36*2,''f'',[1 120000],''nf'',100);',ic)
c_eval('wavE?hmfe.f_units = ''Hz''; wavE?hmfe.f_label = ''f [Hz]''; wavE?hmfe.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)
%c_eval('wavB? = irf_wavelet(gseB?scm.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
%c_eval('wavB?.f_units = ''nT''; wavB?.f_label = ''f [Hz]''; wavB?.p_label = {''log_{10} B^2'',''(nT)^2/Hz''};',ic)

%% Plot detailed electron distribution for selected times
% Initialize figure
[h1,h2] = initialize_combined_plot(10,2,2,0.4,'vertical');
cmap = 'jet';
ic = 4;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:32.10Z');

% Plot
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 1 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1200 700];  
end
if 1 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  ylabel(hca,{'v_{e,\perp}','(km/s)'},'interpreter','tex');
end
if 0 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
end
if 0 % eDist omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.e64.omni.tlim(tint).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [20 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  hca.CLim = [7.3 8.3];
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1par.tlim(tint),gseVe2par.tlim(tint),gseVe3par.tlim(tint),gseVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 1 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E par
  hca = irf_panel('E par');
  try
    set(hca,'ColorOrder',mms_colors('b12'))
    c_eval('irf_plot(hca,{E?par.x,E?par.y,gseE?par},''comp'');',ic)
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    set(hca,'ColorOrder',mms_colors('b12'))
    irf_legend(hca,{'error','E_{||}','E_{||}'},[0.98 0.9],'fontsize',12);
  catch
    set(hca,'ColorOrder',mms_colors('2'))
    c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    set(hca,'ColorOrder',mms_colors('2'))
    irf_legend(hca,{'E_{||}'},[0.98 0.9],'fontsize',12);
  end
  %irf_zoom(hca,'y')
end
if 1 % E
  hca = irf_panel('E hmfe');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?hmfe.x.tlim(tint),dslE?hmfe.y.tlim(tint),dslE?hmfe.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E (DSL)','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end

if 1 % B scm
  hca = irf_panel('B scm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?scm.x.tlim(tint),gseB?scm.y.tlim(tint),gseB?scm.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B (SCM)','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Q
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{Q1.tlim(tint).abs,Q2.tlim(tint).abs,Q3.tlim(tint).abs,Q4.tlim(tint).abs},'comp');    
  hca.YLabel.String = 'Q';
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end

irf_zoom(h1,'x',tintZoom)
irf_zoom(h1([1:4 6:9]),'y')
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);

%% Plot single time particle distributions, 1 sc, 4 projections,
time = irf_time('2015-10-16T10:34:32.83Z','utc>epochtt');
time = irf_time('2015-10-16T11:28:36.40Z','utc>epochtt');
time = irf_time('2015-09-15T11:20:31.74Z','utc>epochtt');
time = irf_time('2015-11-12T07:19:21.30Z','utc>epochtt');
if exist('hmark'); delete(hmark); end
hmark = irf_pl_mark(h1,time.epochUnix','green');
  
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)

% Projection coordinate system
c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
c_eval('hatExB = cross(hatE,hatB);',ic)

y = hatE;
z = hatB;
x = cross(y,z)/norm(cross(y,z));

x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

timeUTC = time.utc;      
isub = 1;

vectors = {hatExB,'ExB'; hatE,'E';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z';v1,'maxvar'};

% Perpendicular plane    
hca = h2(isub); isub = isub + 1; 
xyz = [x;y;z]; vlabels = {'v_{ExB}','v_E','v_B'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
%mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
hca.Title.String = '';
hca.Title.String = timeUTC(1:23);

% B plane 1
hca = h2(isub); isub = isub + 1; 
xyz = [x;z;-y]; vlabels = {'v_{ExB}','v_B','-v_{E}'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
hca.Title.String = '';

% B plane 2
hca = h2(isub); isub = isub + 1;
xyz = [y;z;x]; vlabels = {'v_E','v_B','v_{ExB}'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
hca.Title.String = '';

% Perp plane, but with markings
vectors = {hatExB,'ExB'; hatE,'E';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z';v1,'maxvar'};
hca = h2(isub); isub = isub + 1; 
xyz = [x;y;z]; vlabels = {'v_{ExB}','v_E','v_B'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
hca.Title.String = '';

for ii = 1:4
  colormap(h2(ii),strCMap)
end

%% Plot single time particle distributions, 1 sc, 4 projections,
  
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
c_eval('gradPeslow = gseGradPe.resample(gseVe?);',ic)

% Projection coordinate system

%y = hatE;
%z = hatB;
%x = cross(y,z)/norm(cross(y,z));

x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

c_eval('times = ePDist?.time;',ic)
tind = 886:1:920;953; 911;
tind = 1190:1260;
%tind = 920;
for it = tind
  time = times(it);
  if exist('hmark'); delete(hmark); end
  hmark = irf_pl_mark(h1,time.epochUnix','green');
  
  c_eval('hatGradPe = gradPeslow.resample(time).data/gradPeslow.resample(time).abs.data;',ic)
  c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
  c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
  c_eval('hatExB = cross(hatE,hatB);',ic)
  par = hatB;
  perp1 = hatExB;
  perp2 = cross(par,perp1);
  
  c_eval('hatV = dbcsVe?(it).data/dbcsVe?(it).abs.data;',ic)

  

  timeUTC = time.utc;      
  isub = 1;

  % Perpendicular plane    
  hca = h2(isub); isub = isub + 1; 
  vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B';hatV,'v_e';hatGradPe,'\nabla P_e'};
  xyz = [x;y;z]; vlabels = {'x_{GSE}','y_{GSE}','z_{GSE}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
  hca.Title.String = timeUTC(1:23);

  % B plane 1
  hca = h2(isub); isub = isub + 1; 
  xyz = [x;z;-y]; vlabels = {'x_{GSE}','z_{GSE}','-y_{GSE}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % B plane 2
  hca = h2(isub); isub = isub + 1;
  xyz = [y;z;x]; vlabels = {'y_{GSE}','z_{GSE}','x_{GSE}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % ExB plane, with markings
  vectors = {hatExB,'ExB'; hatE,'E';hatB,'B';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z';hatV,'v_e'};
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{E,\perp}','v_B'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
  hca.Title.String = '';

  for ii = 1:4
    colormap(h2(ii),strCMap)
  end
  cn.print(['e_proj_fix_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','hard_path',eventPath)
end

%%
h = irf_plot(2);
cmap = 'jet';
ic = 4;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');

% Plot
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1
  hca = irf_panel('E spectrogram');
  c_eval('tE = irf_time(wavE?hmfe.t,''epoch>epochtt'');',ic)
  iE = tE.tlim(tintZoom);
  c_eval('plotWavE = wavE?hmfe; plotWavE.t = plotWavE.t(iE); plotWavE.p{1} = plotWavE.p{1}(iE,:);',ic)
  c_eval('irf_spectrogram(hca,plotWavE,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h(:),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);