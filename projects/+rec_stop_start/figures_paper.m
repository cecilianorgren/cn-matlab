tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:09:30.00Z');

%% MMS SWT: Figure: off-time current sheet structure, current
ic = 1;
ic = 1;

comps = ['x','y','z'];
comps = ['y'];

fontsize_leg = 14;
fontsize = 14;

npanels = 1+numel(comps)*2;
nrows = 2;
ncols = 1;
%[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical');
h1 = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
dt_resample = 0.5;
timeline = tint(1):dt_resample:tint(2);
doResample = 1;
isub = 0;
zoomy = [];
matlab_colors = pic_colors('matlab');
j_colors = [mms_colors('123'); matlab_colors(6,:)];
j_colors = [mms_colors('123'); matlab_colors(3,:)];
j_colors = [mms_colors('xyz'); 0 0 0];

plot_colors = [0.8 0.8 0.8; j_colors([4 3],:)];
leg_colors = [0.6 0.6 0.6; j_colors([4 3],:)];
%j_colors = matlab_colors;

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.02 0.9],'fontsize',fontsize_leg);
end 
if 0 % J, Jeav, Jiav, Jav ,curl
  for comp = ['x','y','z']      
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',j_colors)
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize_leg)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',j_colors)
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp)},[0.02 0.99],'fontsize',fontsize_leg);  
  end
end
if 1 % Only  Jav ,curl
  for comp = comps      
    isub = isub + 1;    
    hca = irf_panel(sprintf('J%s',comp));    
    set(hca,'ColorOrder',plot_colors)
    irf_plot(hca,{gseJav.(comp),gseJav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');        
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};    
    set(hca,'ColorOrder',leg_colors)    
    irf_legend(hca,{sprintf('<J_{%s}^{FPI}> (%g ms)',comp,cadence*1e3),sprintf('<J_{%s}^{FPI}> (%g s)',comp,dt_resample),sprintf('J_{%s}^{curl} (%g s)',comp,dt_resample)}',[0.01 0.99],'fontsize',fontsize_leg);  
    %irf_legend(hca,{sprintf('<J_{%s}^{FPI}> (%g ms)',comp,cadence*1e3),sprintf('<J_{%s}^{FPI}> (%g s)',comp,dt_resample),sprintf('J_{%s}^{curl} (%g s)',comp,dt_resample)},[0.02 0.99],'fontsize',fontsize_leg);      
    %hl = findobj(hca,'type','line');    
    %c_eval('hl(?).LineWidth = 0.5;',1:2)
  end  
end
if 1 % Only  Jiav Jeav
  for comp = comps
    isub = isub + 1;    
    hca = irf_panel(sprintf('J%s_ie',comp));            
    set(hca,'ColorOrder',plot_colors)
    irf_plot(hca,{gseJeav.(comp),gseJeav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),},'comp');    
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',leg_colors)     
    
    irf_legend(hca,{sprintf('<J_{e%s}^{FPI}> (%g ms)',comp,cadence*1e3),sprintf('<J_{e%s}^{FPI}> (%g s)',comp,dt_resample)}',[0.01 0.99],'fontsize',fontsize_leg);          
    set(hca,'ColorOrder',leg_colors(3,:))     
    irf_legend(hca,{sprintf('<J_{i%s}^{FPI}> (%g s)',comp,dt_resample)}',[0.01 0.07],'fontsize',fontsize_leg);  
    %irf_legend(hca,{sprintf('<J_{e%s}^{FPI}> (%g ms)',comp,cadence*1e3),sprintf('<J_{e%s}^{FPI}> (%g s)',comp,dt_resample),sprintf('<J_{i%s}^{FPI}> (%g s)',comp,dt_resample)}',[0.01 0.05],'fontsize',fontsize_leg);  
    %irf_legend(hca,{sprintf('<J_{e%s}^{FPI}> (%g ms)',comp,cadence*1e3),sprintf('<J_{e%s}^{FPI}> (%g s)',comp,dt_resample),sprintf('<J_{i%s}^{FPI}> (%g s)',comp,dt_resample)},[0.02 0.99],'fontsize',fontsize_leg);  
    
    %hl = findobj(hca,'type','line');    
    %c_eval('hl(?).LineWidth = 0.5;',1:2)    
  end
end
if 0 % J, mms1-4,curl
  for comp = ['x','y','z']  
    species = 'i';    
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('123b'))
    %irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    eval(sprintf('irf_plot(hca,{gseJ%s1.(comp),gseJ%s2.(comp),gseJ%s3.(comp),gseJ%s4.(comp),gseJcurl.(comp)},''comp'');',species,species,species,species))
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s,%s}^{mms1}',species,comp),sprintf('J_{%s,%s}^{mms2}',species,comp),sprintf('J_{%s,%s}^{mms3}',species,comp),sprintf('J_{%s,%s}^{mms4}',species,comp),'curl'},[0.02 0.99],'fontsize',12);  
  end
end

irf_zoom(h1,'x',tint)
irf_zoom(h1(zoomy),'y')
irf_plot_axis_align 
%irf_pl_mark(h1(1),tint_harris)

h1(1).YLim = [-18 25]*0.99;
h1(2).YLim = [-15 30]*0.99;
h1(3).YLim = [-20 20]*0.99;

%c_eval('h1(?).Position(1) = 0.08;',1:numel(h1))
c_eval('h1(?).Layer = ''top'';',1:numel(h1))
c_eval('h1(?).FontSize = fontsize;',1:numel(h1))
hl = findobj(gcf,'type','line');
%c_eval('hl(?).LineWidth = 0.5;',1:numel(hl))

%annotation('textarrow',[0.33 0.38],[.18 .23],'string',{'<J_{ey}^{FPI}> changes','direction'}','fontsize',fontsize_leg,'horizontalalignment','center');
%annotation('textarrow',[0.36 0.38],[.18 .23],'string',{'<J_{ey}^{FPI}> changes direction','<J_{ey}^{FPI}> remains positive'}','fontsize',fontsize_leg,'horizontalalignment','center');

if 0 % binnin plots of current
isub = 1;
if 0 % Jy,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJav.tlim(tint_harris).y.data;
  yy = gseJcurl.resample(gseJeav).tlim(tint_harris).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 1 % Jy,Jicurl, binning
  hca = h2(isub); isub = isub + 1;
  timeline = gseJeav.time;
  cadence = timeline(2)-timeline(1);
  xx = gseJav.tlim(tint_harris).y.data;
  yy = gseJcurl.resample(timeline).tlim(tint_harris).y.data;
  
  p = polyfit(xx,yy,1);
  %plot(hca,xx,yy,'.')
  jmin = -17;
  jmax = -jmin;
  jstep = 0.5;
  j_edges = jmin:jstep:jmax;
  j_center = (jmin+jstep*0.5):jstep:(jmax-jstep*0.5);
  [N,edges,mid,loc] = histcn([xx yy], j_edges, j_edges);
  
  hca.XLabel.String = '<J_{y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  pcolor(hca,j_center,j_center,(N)')
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hcb.YLabel.String = 'counts';
  hca.CLim = [0 35];
  %plot(hca,[-100 100],p(2)+[-100 100]*p(1))  
  plot(hca,j_edges,j_edges,'color',[0.5 0.5 0.5],'linestyle','--')
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = [jmin jmax];
  hca.YLim = [jmin jmax]; 
  axis(hca,'square')
  drawnow
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  %axis(hca,'equal')
  irf_legend(hca,sprintf('%g ms',cadence*1e3),[0.02 0.98],'fontsize',fontsize_leg)
end

if 1 % Jey,Jiy, binning
  hca = h2(isub); isub = isub + 1;
  timeline = gseJeav.time;
  cadence = timeline(2)-timeline(1);
  xx = gseJeav.tlim(tint_harris).y.data;
  yy = gseJiav.resample(timeline).tlim(tint_harris).y.data;  
  p = polyfit(xx,yy,1);
  %plot(hca,xx,yy,'.')
  jmin = -17;
  jmax = -jmin;
  jstep = 0.5;
  j_edges = jmin:jstep:jmax;
  j_center = (jmin+jstep*0.5):jstep:(jmax-jstep*0.5);
  [N,edges,mid,loc] = histcn([xx yy], j_edges, j_edges);
  
  hca.XLabel.String = '<J_{ey}^{FPI}>';
  hca.YLabel.String = '<J_{iy}^{FPI}>';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  pcolor(hca,j_center,j_center,(N)')
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hcb.YLabel.String = 'counts';
  hca.CLim = [0 35];
  %plot(hca,[-100 100],p(2)+[-100 100]*p(1))  
  plot(hca,j_edges,-j_edges,'color',[0.5 0.5 0.5],'linestyle','--')
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = [jmin jmax];
  hca.YLim = [jmin jmax]; 
  axis(hca,'square')
  drawnow
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  %axis(hca,'equal')   
  irf_legend(hca,sprintf('%g ms',cadence*1e3),[0.02 0.98],'fontsize',fontsize_leg)
  irf_legend(hca,{'data from';'yellow-shaded';'interval'},[0.98 0.08],'fontsize',fontsize_leg,'color',[0 0 0])
end


c_eval('h2(?).FontSize = fontsize;',1:numel(h2))
c_eval('h2(?).Box = ''on'';',1:numel(h2))
end

%% MMS SWT: Figure: off-time current sheet structure, vi vExB
ic = 1;

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

npanels = 6;
%npanels = 1+numel(comps)*2;
nrows = 2;
ncols = 1;
%[h,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical');
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

j_colors = [mms_colors('xyz'); 0 0 0];

plot_colors = [0.8 0.8 0.8; j_colors([4 3],:)];
leg_colors = [0.6 0.6 0.6; j_colors([4 3],:)];

fontsize_leg = 14;

isub = 0;
zoomy = [];
cs = 'gse';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize_leg);
end 
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i\perp,L}','v_{i\perp,M}','v_{i\perp,N}','v_{i||}'}',[1.02 0.9],'fontsize',fontsize_leg);
end
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',fontsize_leg);  
end
if 0 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,1/0.3},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-3:1:2];
  hca.YLabel.Interpreter = 'tex';
end
if 0 % E cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sE?.x,%sE?.y,%sE?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E cs, resamp vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sE?.x.resample(gseVi?),%sE?.y.resample(gseVi?),%sE?.z.resample(gseVi?)},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize_leg);
  irf_zoom(hca,'y')
end

elim_feeps = [8e4 Inf];
if 1 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',plot_colors)
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    c_eval(sprintf('ne = ne?.resample(%sVi?);',cs),ic);    
    %ve.data(abs(ve.data)>5000) = NaN;
    ve.data(abs(ne.data)<0.04,:) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve.resample(%sVi?),%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp%s}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',leg_colors)
    irf_legend(hca,{'v_e','v_i','ExB'}',[1.02 0.9],'fontsize',fontsize_leg);    
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','M','N'];
  comps = ['x','y','z'];
  for icomp = 1:3 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval(sprintf('irf_plot(hca,{%sVe?perp.resample(%sVi?).(comp),%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[1.02 0.9],'fontsize',fontsize_leg);    
    %irf_legend(hca,{'v_e','v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end

if 0 % E+vxB.x/y/z , Vi resample, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve
    isub = isub + 1;    
    hca = irf_panel(['E + VxB ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('evexb = %sE?perp.(comp).resample(%sVi?)+gseVexB?.(comp).resample(%sVi?);',cs,cs,cs),ic)
    c_eval(sprintf('irf_plot(hca,{evexb,%sE?perp.(comp).resample(%sVi?)+gseVixB?.(comp)},''comp'');',cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('E+v_{(%s)}\times B',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'e','i'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
    hca.YLim = [-5 5];
  end
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  eval(sprintf('irf_plot(hca,{%sJcurl.x,%sJcurl.y,%sJcurl.z},''comp'');',cs,cs,cs))
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
if 0 % Vi par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
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

if 0 % EIS Pitch
  isub = isub + 1;
  hca = irf_panel('eis Pitch');  
  irf_spectrogram(hca,eis_pa2.elim([14 4500]*1e3).specrec('pitchangle'));  
  hca.YLim = [0 180];
  %hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.specrec(''energy''),''log'');',2)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',2)
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

if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
  hca.CLim = [-5 -1];
end

if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Te?par/Te?perp},''comp'');',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('xyza'))  
  hca.YLim = [0 5];
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0],'fontsize',fontsize_leg)
  nInd = nInd + 1;
  h(ii).FontSize = 14;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);


hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [0 4]*0.99;

%
hca = irf_panel(['V ExB Vi ','x']); hca.YLim = [-100 250]; hca.YLim = [-150 450]*0.99;
hca = irf_panel(['V ExB Vi ','y']); hca.YLim = [-300 700]*0.99;
hca = irf_panel(['V ExB Vi ','z']); hca.YLim = [-300 400]*0.99;
%hca = irf_panel('E gsm'); hca.YLim = [-15 15]*0.99;

c_eval('h(?).FontSize = 16;',1:numel(h))

% text: 'Moderate beta...'
if 0 % pressure panel included
  %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]-dy);
elseif 0
  annotation('textarrow',[0.3 0.3],[.7 .72],'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]);
  annotation('textarrow',[0.63 0.63],[.7 .72],'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]);  
elseif 0
  %%
  dy = 0.08;
  dx = 0.07; 
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63]-dx,[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63]-dx,[.68 .64]-dy);
end

drawnow
%h(6).YLim(1) = h(7).YLim(2);
%h(2).YLim(1) = 0.001;

drawnow
%hmark = irf_pl_mark(h,time_mark,'k');
%c_eval('hmark(?).LineStyle = '':'';',1:5)
%c_eval('hmark(?).LineWidth = 1;',1:5)

%hat=annotation('textarrow',[0.593 0.593],[.2 .22],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
%hat=annotation('textarrow',[0.593 0.593],[.245 .26],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
if 0 % binning plots of current
isub = 1;
if 0 % ve, vExB scatter
  hca = h2(isub); isub = isub + 1;
  c_eval('xx = gseVe?.tlim(tint_harris).y.data;',ic)
  c_eval('yy = gseVExB?.resample(gseVe?).tlim(tint_harris).y.data;',ic)
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = 'v_{ey}^{FPI}';
  hca.YLabel.String = 'v_{ExB}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  %plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  %plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  %hca.XLim = 40*[-1 1];
  %hca.YLim = 40*[-1 1];  
end
if 1 % ve, vExB, binning
  hca = h2(isub); isub = isub + 1;
  timeline = gseJeav.time;
  cadence = timeline(2)-timeline(1);
  c_eval('xx = gseVe?.tlim(tint_harris).y.data;',ic)
  c_eval('yy = gseVExB?.resample(gseVe?).tlim(tint_harris).y.data;',ic)
  
  p = polyfit(xx,yy,1);
  %plot(hca,xx,yy,'.')
  jmin = -1000;
  jmax = -jmin;
  jstep = 20;
  j_edges = jmin:jstep:jmax;
  j_center = (jmin+jstep*0.5):jstep:(jmax-jstep*0.5);
  [N,edges,mid,loc] = histcn([xx yy], j_edges, j_edges);
  
  hca.XLabel.String = 'v_{ey}^{FPI}';
  hca.YLabel.String = 'v_{ExB}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  pcolor(hca,j_center,j_center,log10(N)')
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hcb.YLabel.String = 'counts';
  %hca.CLim = [0 35];
  %plot(hca,[-100 100],p(2)+[-100 100]*p(1))  
  plot(hca,j_edges,j_edges,'color',[0.5 0.5 0.5],'linestyle','--')
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = [jmin jmax];
  hca.YLim = [jmin jmax]; 
  axis(hca,'square')
  drawnow
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  %axis(hca,'equal')
  irf_legend(hca,sprintf('%g ms',cadence*1e3),[0.02 0.98],'fontsize',fontsize_leg)
end
if 1 % vi, vExB, binning
  hca = h2(isub); isub = isub + 1;
  timeline = gseJeav.time;
  cadence = timeline(2)-timeline(1);
  c_eval('xx = gseVi?.tlim(tint_harris).y.data;',ic)
  c_eval('yy = gseVExB?.resample(gseVi?).tlim(tint_harris).y.data;',ic)
  
  p = polyfit(xx,yy,1);
  %plot(hca,xx,yy,'.')
  jmin = -1000;
  jmax = -jmin;
  jstep = 20;
  j_edges = jmin:jstep:jmax;
  j_center = (jmin+jstep*0.5):jstep:(jmax-jstep*0.5);
  [N,edges,mid,loc] = histcn([xx yy], j_edges, j_edges);
  
  hca.XLabel.String = 'v_{iy}^{FPI}';
  hca.YLabel.String = 'v_{ExB}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  pcolor(hca,j_center,j_center,log10(N)')
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hcb.YLabel.String = 'counts';
  %hca.CLim = [0 35];
  %plot(hca,[-100 100],p(2)+[-100 100]*p(1))  
  plot(hca,j_edges,j_edges,'color',[0.5 0.5 0.5],'linestyle','--')
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = [jmin jmax];
  hca.YLim = [jmin jmax]; 
  axis(hca,'square')
  drawnow
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  %axis(hca,'equal')
  irf_legend(hca,sprintf('%g ms',cadence*1e3),[0.02 0.98],'fontsize',fontsize_leg)
end


c_eval('h2(?).FontSize = fontsize;',1:numel(h2))
c_eval('h2(?).Box = ''on'';',1:numel(h2))
end
