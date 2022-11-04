tint_b = irf.tint('2017-07-25T22:09:00.82Z/2017-07-25T22:11:31.82Z');

tint_c = irf.tint('2017-07-25T22:09:00.82Z/2017-07-25T22:11:00.00Z');
tint_d = irf.tint('2017-07-25T22:09:30.82Z/2017-07-25T22:11:00.00Z');

tint_action = irf.tint('2017-07-25T22:09:40.00Z/2017-07-25T22:10:50.00Z');
tint_cavity = irf.tint('2017-07-25T22:09:46.00Z/2017-07-25T22:09:56.00Z');

tint_ps = irf.tint('2017-07-25T22:09:45.00Z/2017-07-25T22:10:06.00Z');
tint_figure = irf.tint('2017-07-25T22:09:25.00Z/2017-07-25T22:10:40.00Z');
tint_figure = irf.tint('2017-07-25T22:09:15.00Z/2017-07-25T22:10:40.00Z');

%% Prepare data

% Integrate vExB and vi over the time interval of interest
c_eval('intVex? = irf_integrate(gseVe?perp.x,tint_ps(1));',ic)
c_eval('intVix? = irf_integrate(gseVi?perp.x,tint_ps(1));',ic)
c_eval('intVExBx? = irf_integrate(gseVExB?.x,tint_ps(1));',ic)

% What is the different over the designated interval
c_eval('dL_vex? = diff(intVex?.resample(tint_ps).data);',ic)
c_eval('dL_vix? = diff(intVix?.resample(tint_ps).data);',ic)
c_eval('dL_vExBx? = diff(intVExBx?.resample(tint_ps).data);',ic)

% LMN

% Integrate vExB and vi over the time interval of interest
c_eval('lmnintVex? = irf_integrate(lmnVe?perp.x,tint_ps(1));',ic)
c_eval('lmnintVix? = irf_integrate(lmnVi?perp.x,tint_ps(1));',ic)
c_eval('lmnintVExBx? = irf_integrate(lmnVExB?.x,tint_ps(1));',ic)

% What is the different over the designated interval
c_eval('dL_vex? = diff(lmnintVex?.resample(tint_ps).data);',ic)
c_eval('dL_vix? = diff(lmnintVix?.resample(tint_ps).data);',ic)
c_eval('dL_vExBx? = diff(lmnintVExBx?.resample(tint_ps).data);',ic)
%% Prepare reduced ion distribution for shorter time interval
vint = [-Inf Inf];
elim_fig = [800 50000];
c_eval('par = dmpaB?.norm.resample(iPDist?);',ic)
L = R(1,:);
c_eval('if1DL? = iPDist?.tlim(tint_figure+[-5 5]).elim(elim_fig).reduce(''1D'',L,''vint'',vint);',ic)
c_eval('if1Dpar? = iPDist?.tlim(tint_figure+[-5 5]).elim(elim_fig).reduce(''1D'',par,''vint'',vint);',ic)
c_eval('elim_effective = iPDist?.elim(elim_fig).ancillary.energy(1,1) - iPDist?.elim(elim_fig).ancillary.delta_energy_minus(1,1);',ic)

%% Figure: 
ic = 1;

npanels = 7;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));

isub = 0;
zoomy = [];

if 0 % B gs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z,gseB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end 
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % vExB gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmVExB?.x,gsmVExB?.y,gsmVExB?.z},''comp'');',ic)  
  c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 1 % vExB gse, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.resample(gseVi?).x,gseVExB?.resample(gseVi?).y,gseVExB?.resample(gseVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 0 % J 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x,gseJ?.y,gseJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 0 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?perp.x,gseVi?perp.y,gseVi?perp.z,gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
if 1 % Ve perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
 
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 0.5;
  %c_eval('irf_plot(hca,{gseVe?perp.x.filt(0,ffilt,[],5),gseVe?perp.y.filt(0,ffilt,[],5),gseVe?perp.z.filt(0,ffilt,[],5),gseVe?par.filt(0,ffilt,[],5)},''comp'');',ic)  
  c_eval('irf_plot(hca,{gseVe?perp.x.resample(iPDist?),gseVe?perp.y.resample(iPDist?),gseVe?perp.z.resample(iPDist?),gseVe?par.resample(iPDist?)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 1 % vExB gse, Vi  gse, compare x
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB vi x');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.resample(gseVi?).x,gseVi?perp.x,gseVe?perp.x.resample(gseVi?)},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{x}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{ExB}','v_{i,\perp}','v_{e,\perp}'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 1 % time integrated vExB gse, Vi  gse, compare x
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('int V ExB vi x');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{intVExBx?,intVix?,intVex?},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'L = \int v_{x} dt','(km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{ExB}','v_{i,\perp}','v_{e,\perp}'},[0.98 0.9],'fontsize',12);
  
  c_eval('dL_vex = dL_vex?;',ic)
  c_eval('dL_vix = dL_vix?;',ic)
  c_eval('dL_vExBx = dL_vExBx?;',ic)
  irf_legend(hca,{sprintf('L^{v_{ExB}} = %.0f km',dL_vExBx),sprintf('L^{v_i} = %.0f km',dL_vix),sprintf('L^{v_e} = %.0f km',dL_vex)}',[0.1 0.9],'fontsize',12);
  
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint_figure)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')  
hca = irf_panel('V ExB vi x'); hca.YLim = [-200 2000];
hca = irf_panel('int V ExB vi x'); hca.YLim(1) = -5e3; hca.YLim(2) = 5.5e4;
hca = irf_panel('Ve perp par'); hca.YLim = [-3e3 3e3];

irf_plot_axis_align
irf_pl_mark(h,tint_ps)
h(1).Title.String = irf_ssub('MMS ?',ic);

%% Figure: compact
ic = 1;

npanels = 7;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));

fontsize = 14;
isub = 0;
zoomy = [];

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z,gseB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.05 0.7],'fontsize',fontsize);  
end 
if 1 % i psd x
  isub = isub + 1;
  hca = irf_panel('iLine x');
  c_eval('if1D = if1Dx?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    %irf_plot(hca,gseVi1.x,'k')
    irf_plot(hca,lmnVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i psd L
  isub = isub + 1;
  hca = irf_panel('iLine L');
  c_eval('if1D = if1DL?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    %irf_plot(hca,gseVi1.x,'k')
    irf_plot(hca,lmnVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i psd x f(vx)*vz
  isub = isub + 1;
  hca = irf_panel('iLine x flux');
  c_eval('if1D = if1D?M;',ic)
  c_eval('if1D = if1Dx?; if1D.data = if1D.data.*if1D.depend{1}*1e3; if1D.units = ''1/m^3'';',ic)  % right units?
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 1 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.x,'k')
    %irf_plot(hca,lmnVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
if 1 % i psd x f(vx)*vz
  isub = isub + 1;
  hca = irf_panel('iLine x flux L');
  c_eval('if1D = if1D?M;',ic)
  c_eval('if1D = if1DL?; if1D.data = if1D.data.*if1D.depend{1}*1e3; if1D.units = ''1/m^3'';',ic)  % right units?
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 1 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.x,'k')
    %irf_plot(hca,lmnVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % vExB gse, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.resample(gseVi?).x,gseVExB?.resample(gseVi?).y,gseVExB?.resample(gseVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?perp.x,gseVi?perp.y,gseVi?perp.z,gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
 
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 0.5;
  %c_eval('irf_plot(hca,{gseVe?perp.x.filt(0,ffilt,[],5),gseVe?perp.y.filt(0,ffilt,[],5),gseVe?perp.z.filt(0,ffilt,[],5),gseVe?par.filt(0,ffilt,[],5)},''comp'');',ic)  
  c_eval('irf_plot(hca,{gseVe?perp.x.resample(iPDist?),gseVe?perp.y.resample(iPDist?),gseVe?perp.z.resample(iPDist?),gseVe?par.resample(iPDist?)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 1 % vExB gse, Vi  gse, compare x
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB vi x');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.resample(gseVi?).x,gseVi?perp.x,gseVe?perp.x.resample(gseVi?)},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,x}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{ExB}','v_{i,\perp}','v_{e,\perp}'},[0.05 0.8],'fontsize',fontsize);
  hca.Children = hca.Children(end:-1:1);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 1 % time integrated vExB gse, Vi  gse, compare x
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('int V ExB vi x');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{intVExBx?*1e-3,intVix?*1e-3,intVex?*1e-3},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'L = \int v_{\perp,x} dt','(10^3 km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{ExB}','v_{i,\perp}','v_{e,\perp}'},[0.98 0.9],'fontsize',12);
  
  c_eval('dL_vex = dL_vex?;',ic)
  c_eval('dL_vix = dL_vix?;',ic)
  c_eval('dL_vExBx = dL_vExBx?;',ic)
  irf_legend(hca,{sprintf('L^{ExB} = %.0f km',100*round(dL_vExBx/100)),sprintf('L^{i} = %.0f km',100*round(dL_vix/100)),sprintf('L^{e} = %.0f km',100*round(dL_vex/100))}',[0.05 0.8],'fontsize',fontsize);
  
end

if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.elim(elim_fig).omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
irf_legend(h(1),legends{1},[0.01 0.98],'color',[0 0 0],'fontsize',16)
irf_legend(h(2),legends{2},[0.01 0.98],'color',[0 0 0],'fontsize',16)
irf_legend(h(3),legends{3},[0.01 0.98],'color',[0 0 0],'fontsize',16)
for ii = 1:npanels  
  nInd = nInd + 1;
  h(ii).FontSize = fontsize;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint_figure)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')  
hca = irf_panel('V ExB vi x'); hca.YLim = [-300 2900];
hca = irf_panel('int V ExB vi x'); hca.YLim(1) = -5e3*1e-3; hca.YLim(2) = 5.5e4*1e-3;
%hca = irf_panel('Ve perp par'); hca.YLim = [-3e3 3e3];

irf_plot_axis_align
irf_pl_mark(h,tint_ps)
h(1).Title.String = irf_ssub('MMS ?',ic);
%c_eval('h(?).FontSize = 14;',1:npanels)

%% Figure: compact, LMN
ic = 1;

npanels = 4;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));

set(gcf,'position',[ 611   199   670   498]); drawnow
fontsize = 14;
isub = 0;
zoomy = [];

if 1 % B LMN
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B lmn');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{lmnB?.x,lmnB?.y,lmnB?.z,lmnB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_L','B_M','B_N','|B|'},[0.05 0.7],'fontsize',fontsize);  
  irf_legend(hca,{'B_L','B_M','B_N','|B|'}',[1.02 0.9],'fontsize',fontsize);  
end
if 0 % i psd par, initial field aligned beam slightly more visible in this format
  isub = isub + 1;
  hca = irf_panel('f(vpar)');
  c_eval('if1D = if1Dpar?;',ic)
  if1D.data(if1D.data==0) = NaN;
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));    
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    %irf_plot(hca,gseVi1.x,'k')
    irf_plot(hca,lmnVi1par,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{i||}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,sprintf('E_i > %.0f eV',round(elim_effective)),[0.98 0.08],'k')
end
if 1 % i psd L
  isub = isub + 1;
  hca = irf_panel('f(vL)');
  c_eval('if1D = if1DL?;',ic)
  if1D.data(if1D.data==0) = NaN;
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));    
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    %irf_plot(hca,gseVi1.x,'k')
    irf_plot(hca,lmnVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iL}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,sprintf('E_i > %.0f eV',round(elim_effective)),[0.98 0.08],'k')
end
if 0 % i psd L f(vL)*vL
  isub = isub + 1;
  hca = irf_panel('iLine x flux L');
  c_eval('if1D = if1D?M;',ic)
  c_eval('if1D = if1DL?; if1D.data = if1D.data.*if1D.depend{1}*1e3; if1D.units = ''1/m^3'';',ic)  % right units?
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 1 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.x,'k')
    %irf_plot(hca,lmnVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
end
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % vExB gse, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.resample(gseVi?).x,gseVExB?.resample(gseVi?).y,gseVExB?.resample(gseVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?perp.x,gseVi?perp.y,gseVi?perp.z,gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
 
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 0.5;
  %c_eval('irf_plot(hca,{gseVe?perp.x.filt(0,ffilt,[],5),gseVe?perp.y.filt(0,ffilt,[],5),gseVe?perp.z.filt(0,ffilt,[],5),gseVe?par.filt(0,ffilt,[],5)},''comp'');',ic)  
  c_eval('irf_plot(hca,{gseVe?perp.x.resample(iPDist?),gseVe?perp.y.resample(iPDist?),gseVe?perp.z.resample(iPDist?),gseVe?par.resample(iPDist?)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 1 % vExB LMN, Vi  gse, compare x
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB vi x');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{lmnVExB?.resample(lmnVi?).x,lmnVi?perp.x,lmnVe?perp.x.resample(lmnVi?)},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{{\perp}L}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{ExB}','v_{i,\perp}','v_{e,\perp}'},[0.05 0.8],'fontsize',fontsize);
  hca.Children = hca.Children(end:-1:1);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 1 % time integrated vExB gse, Vi  gse, compare x
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('int V ExB vi x');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{intVExBx?*1e-3,intVix?*1e-3,intVex?*1e-3},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'L = \int v_{\perp,x} dt','(10^3 km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{ExB}','v_{i,\perp}','v_{e,\perp}'},[0.98 0.9],'fontsize',12);
  
  c_eval('dL_vex = dL_vex?;',ic)
  c_eval('dL_vix = dL_vix?;',ic)
  c_eval('dL_vExBx = dL_vExBx?;',ic)
  irf_legend(hca,{sprintf('L^{ExB} = %.0f km',100*round(dL_vExBx/100)),sprintf('L^{i} = %.0f km',100*round(dL_vix/100)),sprintf('L^{e} = %.0f km',100*round(dL_vex/100))}',[0.05 0.8],'fontsize',fontsize);
  
end

if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.elim(elim_fig).omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
irf_legend(h(1),legends{1},[0.01 0.98],'color',[0 0 0],'fontsize',16)
irf_legend(h(2),legends{2},[0.01 0.98],'color',[0 0 0],'fontsize',16)
irf_legend(h(3),legends{3},[0.01 0.98],'color',[0 0 0],'fontsize',16)
irf_legend(h(4),legends{4},[0.01 0.98],'color',[0 0 0],'fontsize',16)
for ii = 1:npanels  
  nInd = nInd + 1;
  h(ii).FontSize = fontsize;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint_figure)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')  
hca = irf_panel('V ExB vi x'); hca.YLim = [-300 2900];
hca = irf_panel('int V ExB vi x'); hca.YLim(1) = -5e3*1e-3; hca.YLim(2) = 5.5e4*1e-3;
%hca = irf_panel('Ve perp par'); hca.YLim = [-3e3 3e3];
hca = irf_panel('f(vL)'); hca.CLim = [-5 -0.8]; colormap(hca,pic_colors('pasteljet'))
%hca = irf_panel('f(vpar)'); hca.CLim = [-5 -0.8]; colormap(hca,pic_colors('pasteljet'))

irf_plot_axis_align
irf_pl_mark(h,tint_ps)
%h(1).Title.String = irf_ssub('MMS ?',ic);
%c_eval('h(?).FontSize = 14;',1:npanels)
annotation('textarrow',[0.27 0.27],[.74 .72],'string',{'field-aligned','ions appearing'}','fontsize',12,'horizontalalignment','center');
annotation('textarrow',[0.37 0.40],[.86 .86],'string',{'magnetic topology','changing'}','fontsize',12,'horizontalalignment','center');
annotation('textarrow',[0.575 0.575],[.57 .63],'string',{'cold ions'}','fontsize',12,'horizontalalignment','center');
%% Figure: V ExB
ic = 2;

Etop_fpi = iPDist2.ancillary.energy(1,end)+iPDist2.ancillary.delta_energy_plus(1,end);

npanels = 6;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));

isub = 0;
zoomy = [];

if 0 % B gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'PDi?,',...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_{dyn,i}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 0 % vExB gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.x,gsmVExB?.y,gsmVExB?.z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 1 % vExB gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x,gseVExB?.y,gseVExB?.z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVExB?.resample(gsmVi?).y,gsmVExB?.resample(gsmVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi x');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVi?perp.x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,x}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi y');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).y,gsmVi?perp.y},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,y}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi z');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).z,gsmVi?perp.z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,z}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 0 % J 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x,gseJ?.y,gseJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 0 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?perp.x,gseVi?perp.y,gseVi?perp.z,gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve gsm
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve filt
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 3;
  c_eval('irf_plot(hca,{gsmVe?.filt(0,3,[],3).x.tlim(tint),gsmVe?.filt(0,3,[],3).y.tlim(tint),gsmVe?.filt(0,3,[],3).z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve  gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x,gsmVe?.y,gsmVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x,gseVe?perp.y,gseVe?perp.z},''comp'');',ic)  
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
if 0 % Ve x B
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x,gseVExB?.y,gseVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'-v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',ic)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
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
if 0 % Ti par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTi?.xx.tlim(tint),(facTi?.yy+facTi?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_i','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{i,||}','T_{i,\perp}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end

if 0 % i psd x,y,z, 3 panels
  for comp = ['x','y','z']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 1 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
h(1).Title.String = irf_ssub('MMS ?',ic);

%% Figure: momentum equation / ohms law
ic = 2;

npanels = 7;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));

isub = 0;
zoomy = [];

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z,gseB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end 
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, 4sc average
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_plot(hca,{PBav,Piav,Peav,PDiav,Piav + Peav.resample(Piav) + PBav.resample(Piav)},'comp')
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i}','P_{e}','P_{i,dyn,x}','P_{B}+P_i+P_e'},[0.98 0.9],'fontsize',12);  
end
if 0 % Pressures, single sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'PDi?,',...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_{dyn,i}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 1 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');  
  hca.YLabel.String = {'J^{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 0 % JxB, 3 panels: x,y,z
  for comp = ['x','y','z']
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    scale = 1e0*1e0; units = 'T Am^{-2}';
    scale = 1e9*1e9; units = 'nT nAm^{-2}';
    
    hca = irf_panel(['JxB ', comp]);
    set(hca,'ColorOrder',pic_colors('matlab'))   
    irf_plot(hca,{gseJxB.(comp)*scale,-gseDivPb.(comp)*scale,-gseCurvB.(comp)*scale},'comp');  
    %irf_plot(hca,{gseJxB.(comp)*scale,-gseDivPb.(comp)*scale,-gseCurvB.(comp)*scale,-gseDivPb.(comp)*scale-gseCurvB.(comp)*scale},'comp');  
    hca.YLabel.String = {['JxB_' comp],['(' units ')']};
    set(hca,'ColorOrder',pic_colors('matlab'))
    irf_legend(hca,{'JxB','-\nabla B^2/2\mu_0','-\mu_0^{-1}B\cdot\nabla B'}',[1.01 0.9],'fontsize',12);  
  end
end
if 0 % electron mom eq. 3 panels: x,y,z
  for comp = ['x','y','z']
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    hca = irf_panel(['ele eom ', comp]);
    set(hca,'ColorOrder',mms_colors('xyza'))
    to_plot_gseVexBav = gseVexBav;
    to_plot_gseVexBav.data(abs(to_plot_gseVexBav.data)>100) = NaN;
    irf_plot(hca,{gseEav.(comp),-1*to_plot_gseVexBav.(comp),-1*gseGradPene.resample(gseVi1).(comp),-1*(gseEav.(comp).resample(to_plot_gseVexBav)+to_plot_gseVexBav.(comp))},'comp');  
    hca.YLabel.String = {['E_' comp],'(mV/m)'};
    set(hca,'ColorOrder',mms_colors('xyza'))
    irf_legend(hca,{'E','-v_exB','-divPe/ne','-(E+v_exB)'}',[1.02 0.9],'fontsize',12);  
  end
end
if 1 % ohm's law, 3 panels: x,y,z
  for comp = ['x','y','z']
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    hca = irf_panel(['ele eom ', comp]);
    set(hca,'ColorOrder',pic_colors('matlab'))
    to_plot_gseVixBav = gseVixBav;
    to_plot_gseVixBav.data(abs(to_plot_gseVixBav.data)>100) = NaN;
    irf_plot(hca,{gseEav.(comp),to_plot_gseVixBav.(comp),-1*gseGradPene.(comp).resample(gseVi1),gseJxBne_mVm.(comp),+1*(gseEav.(comp).resample(to_plot_gseVixBav)+to_plot_gseVixBav.(comp))},'comp');  
    hca.YLabel.String = {['E_' comp],'(mV/m)'};
    set(hca,'ColorOrder',pic_colors('matlab'))
    irf_legend(hca,{'E','v_ixB','-divPe/ne','JxB/ne','(E+v_ixB)'}',[1.01 0.9],'fontsize',12);  
  end
end
if 0 % vExB gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmVExB?.x,gsmVExB?.y,gsmVExB?.z},''comp'');',ic)  
  c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVExB?.resample(gsmVi?).y,gsmVExB?.resample(gsmVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi x');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVi?perp.x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,x}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi y');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).y,gsmVi?perp.y},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,y}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi z');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).z,gsmVi?perp.z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,z}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % J 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x,gseJ?.y,gseJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 0 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?perp.x,gseVi?perp.y,gseVi?perp.z,gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
if 1 % Ve gse
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve gsm
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve filt
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 3;
  c_eval('irf_plot(hca,{gsmVe?.filt(0,3,[],3).x.tlim(tint),gsmVe?.filt(0,3,[],3).y.tlim(tint),gsmVe?.filt(0,3,[],3).z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve  gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x,gsmVe?.y,gsmVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x,gseVe?perp.y,gseVe?perp.z},''comp'');',ic)  
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
if 0 % Ve x B
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x,gseVExB?.y,gseVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'-v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x,gseVe?perp.y,gseVe?perp.z,gseVe?par},''comp'');',ic)  
  hca.YLabel.String = {['v_{e} mms' num2str(ic)],'(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\perp x','\perp y','\perp z','||'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp par av
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{gseVeperpav.x,gseVeperpav.y,gseVeperpav.z,gseVeparav},'comp');
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\perp x','\perp y','\perp z','||'},[0.98 0.9],'fontsize',12);
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

irf_zoom(h,'x',tint_action)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%tmark_eis = tint(1):20:tint(2);

%% Figure: V ExB
ic = 2;

Etop_fpi = iPDist2.ancillary.energy(1,end)+iPDist2.ancillary.delta_energy_plus(1,end);

npanels = 8;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));

isub = 0;
zoomy = [];

if 1 % B gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 0 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'PDi?,',...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_{dyn,i}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 0 % vExB gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmVExB?.x,gsmVExB?.y,gsmVExB?.z},''comp'');',ic)  
  c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVExB?.resample(gsmVi?).y,gsmVExB?.resample(gsmVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi x');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVi?perp.x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,x}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi y');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).y,gsmVi?perp.y},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,y}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi z');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).z,gsmVi?perp.z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,z}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 0 % J 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x,gseJ?.y,gseJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 0 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?perp.x,gseVi?perp.y,gseVi?perp.z,gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve gsm
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve filt
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 3;
  c_eval('irf_plot(hca,{gsmVe?.filt(0,3,[],3).x.tlim(tint),gsmVe?.filt(0,3,[],3).y.tlim(tint),gsmVe?.filt(0,3,[],3).z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve  gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x,gsmVe?.y,gsmVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x,gseVe?perp.y,gseVe?perp.z},''comp'');',ic)  
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
if 0 % Ve x B
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x,gseVExB?.y,gseVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'-v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',ic)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
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
if 0 % Ti par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTi?.xx.tlim(tint),(facTi?.yy+facTi?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_i','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{i,||}','T_{i,\perp}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end

if 1 % i psd x,y,z, 3 panels
  for comp = ['x','y','z']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 1 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
h(1).Title.String = irf_ssub('MMS ?',ic);

%% Compare to simulations to see if idea is reasonable
%no02m = PIC('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/fields.h5');

twpe = [3000 22000];
xlim = 80 + 0.2*[-1 1];
zlim = 0.0 + 0.2*[-1 1];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);

t = pic.twci;
Bz = squeeze(mean(mean(pic.Bz,1),2));
Ey = squeeze(mean(mean(pic.Ey,1),2));
vx_ExB = squeeze(mean(mean(pic.vExBx,1),2));
vx_i = squeeze(mean(mean(pic.viperpx,1),2));
vx_e = squeeze(mean(mean(pic.veperpx,1),2));

% Integrate speeds to get some length
L_ExB = cumtrapz(t,vx_ExB);
L_i = cumtrapz(t,vx_i);
L_e = cumtrapz(t,vx_e);

% Fake reconnection rate
R = 0.1;

% Find DF
Blim = 0.5;
iDF = find(abs(Bz)>Blim*max(abs(Bz)),1,'first');
tDF = t(iDF);
inL_ExB = L_ExB(iDF)*R;
inL_vi = L_i(iDF)*R;
inL_ve = L_e(iDF)*R;

% Plot
h = setup_subplots(5,1); isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,t,Bz,tDF,Bz(iDF),'*')
hca.XLabel.String = 't\omega_{ci}';
hca.YLabel.String = 'B_z';

hca = h(isub); isub = isub + 1;
plot(hca,t,Ey)
hca.XLabel.String = 't\omega_{ci}';
hca.YLabel.String = 'E_y';

hca = h(isub); isub = isub + 1;
plot(hca,t,vx_ExB,t,vx_i,t,vx_e)
hca.XLabel.String = 't\omega_{ci}';
hca.YLabel.String = 'v';
legend(hca,{'v_{ExB}','v_i','v_e'},'box','off')
hca.YLim = 1.2*[-1 1];

hca = h(isub); isub = isub + 1;
plot(hca,t,L_ExB,t,L_i,t,L_e)
hca.XLabel.String = 't\omega_{ci}';
hca.YLabel.String = 'L';
legend(hca,{'L_{ExB}','L_i','L_e'},'box','off')


hca = h(isub); isub = isub + 1;
plot(hca,t,L_ExB*R,t,L_i*R,t,L_e*R)
hca.XLabel.String = 't\omega_{ci}';
hca.YLabel.String = 'L*R';
hold(hca,'on')
plot(hca,tDF,inL_ExB,'*',tDF,inL_vi,'*',tDF,inL_ve,'*')
hold(hca,'off')
legend(hca,{'L_{ExB}','L_i','L_e',sprintf('L_{ExB} = %.2f',inL_ExB),sprintf('L_{i} = %.2f',inL_vi),sprintf('L_{e} = %.2f',inL_ve)},...
  'box','off','location','best')

c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on''; h(?).Layer = ''top''; h(?).GridAlpha = 0.1;',1:numel(h))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).FontSize = 14;',1:numel(h))

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

compact_panels(h,0.0)

h(1).Title.String = sprintf('x = [%.1f,%.1f], z = [%.1f,%.1f]',xlim(1),xlim(2),zlim(1),zlim(2));










%% Compare to simulations to see if idea is reasonable, try to show several locations at once
%no02m = PIC('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/fields.h5');

twpe = [0000 25000];
xlim = [65 100];
zlim = 0.0 + 0.2*[-1 1];
pic = no02m.twpelim(twpe).xlim(xlim).zlim(zlim);

t = pic.twci;
Ay = squeeze(mean(pic.A,2));
Bz = squeeze(mean(pic.Bz,2));
Ey = squeeze(mean(pic.Ey,2));
ni = squeeze(mean(pic.ni,2));
vx_ExB = squeeze(mean(pic.vExBx,2));
vx_i = squeeze(mean(pic.viperpx,2));
vx_e = squeeze(mean(pic.veperpx,2));

% Integrate speeds to get some length
L_ExB = cumtrapz(t,vx_ExB,2); % Too noisy
L_i = cumtrapz(t,vx_i,2);
L_e = cumtrapz(t,vx_e,2);
%
%% Fake reconnection rate
R = 0.2;

% Find DF
% Blim = 0.5;
% iDF = find(abs(Bz)>Blim*max(abs(Bz)),1,'first');
% tDF = t(iDF);
% inL_ExB = L_ExB(iDF)*R;
% inL_vi = L_i(iDF)*R;
% inL_ve = L_e(iDF)*R;

% Plot
ALim = [8 12];
LinLim = [0 3]*0.99;

figure(72)
h = setup_subplots(7,1); isub = 1;

if 0 % ni
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,t,ni')
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 't\omega_{ci}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'n_i';
end
if 1 % Ay
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,t,Ay')
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 't\omega_{ci}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'A_y';
  hca.CLim = ALim;
end
if 1 % Bz
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,t,Bz')
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 't\omega_{ci}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'B_z';
end
if 0 % Ey
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,t,Ey')
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 't\omega_{ci}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_y';
  hca.CLim = [0 0.5];
end
if 1 % Vi
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,t,vx_i')
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 't\omega_{ci}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ix\perp}';
end
if 1 % Ve
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,t,vx_e')
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 't\omega_{ci}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'v_{ex\perp}';
end
if 0 % L_vi
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,t,L_i')
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 't\omega_{ci}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'L_{ix\perp}';
end
if 0 % L_ve
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,t,L_e')
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 't\omega_{ci}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'L_{ex\perp}';
end
if 1 % L_vi*R
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,t,abs(R*L_i)')
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 't\omega_{ci}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = '|R*L_{ix\perp}|';
  hca.CLim = LinLim;
end
if 1 % L_ve*R
  hca = h(isub); isub = isub + 1;
  pcolor(hca,pic.xi,t,abs(R*L_e)')
  shading(hca,'flat')
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 't\omega_{ci}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = '|R*L_{ex\perp}|';
  hca.CLim = LinLim;
end
if 1 % L_ve*R
  hca = h(isub); isub = isub + 1;
  %[T,X] = ndgrid(t,pic.xi);
  toplot = abs(R*L_e);
  %toplot(abs(Bz)<0.1) = NaN;
  %toplot(X-100<-1*T) = NaN;
  %toplot = ((X-100)./T)'+1;
  %pcolor(hca,pic.xi,t,toplot')
  %shading(hca,'flat')
  contourf(hca,pic.xi,t,toplot',0:0.2:4)
  hca.XLabel.String = 'x/d_i';
  hca.YLabel.String = 't\omega_{ci}';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = '|R*L_{ex\perp}|';
  hca.CLim = LinLim;
end

for ip = [1 3:numel(h)] %1:numel(h) % Add Bx contours
  hca = h(ip);
  hold(hca,'on')
  clim = hca.CLim;
  contour(hca,pic.xi,t,Bz',[-1:0.2:-0.2 0.2:0.2:1],'color','k','linewidth',1)
  hca.CLim = clim;
  hold(hca,'off')
end
for ip = 2
  hca = h(ip);
  hold(hca,'on')
  clim = hca.CLim;
  contour(hca,pic.xi,t,abs(R*L_e)',0:0.2:4,'color','k','linewidth',1)
  hca.CLim = clim;
  hold(hca,'off')
end
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on''; h(?).Layer = ''top''; h(?).GridAlpha = 0.1;',1:numel(h))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).FontSize = 14;',1:numel(h))

hl = findobj(gcf,'type','line');
%c_eval('hl(?).LineWidth = 1;',1:numel(hl))

compact_panels(h,0.0)
hlinks = linkprop(h,{'XLim','YLim'});

h(1).Title.String = sprintf('z = [%.1f,%.1f]',zlim(1),zlim(2));



%%


figure(73)
if 1 % Ay   
  xx = no02m.twpelim(twpe(1)).x_xline;  
  no02m.twpelim(twpe).xlim(xx+[-0.5 0.5]).zlim([-5 5]).plot_timemap('tz',{'A','ni'}','cmap',{pic_colors('candy4'),pic_colors('candy4')},'A',1,'clim',{ALim,[0 1]})
  %pcolor(hca,pic.xi,t,Ay')
  %shading(hca,'flat')
  %hca.XLabel.String = 'x/d_i';
  %hca.YLabel.String = 't\omega_{ci}';
  %hcb = colorbar('peer',hca);
  %hcb.YLabel.String = 'A_y';
end


