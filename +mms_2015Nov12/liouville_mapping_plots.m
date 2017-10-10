%% Observed fields, for comparison
units = irf_units;
% Observed data
ic = 1;
toffset1 = 0; toffset2 = -0.03; toffset3 = 0.1; toffset4 = 0.1;
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
c_eval('tintObs = tintObs+-toffset?;',ic) % tinObs(2) correspond
CS_normal_velocity = 70; % km/s
%CS_normal_velocity = 70; % km/s

tintObs = tintObs + 1*[-1 1];

c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsCurvB = mvaCurvB.resample(obsB).tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs); obsE = obsE.resample(obsB);'...
'obsE = mvaE?_new.tlim(tintObs); obsE = obsE.resample(obsB);'... %'obsE = mvaE?_new.tlim(tintObs) - mvaEht.resample(mvaE?_new).tlim(tintObs); obsE = obsE.resample(obsB);'...
'obsEpar = mvaE?par.tlim(tintObs); obsEpar = obsEpar.resample(obsB);'...
'obsEperp = mvaE?perp.tlim(tintObs); obsEperp = obsEperp.resample(obsB);'...
'obsVepar = mvaVe?par.tlim(tintObs); obsVepar = obsVepar.resample(obsB);'...
'obsVeperp = mvaVe?perp.tlim(tintObs); obsVeperp = obsVeperp.resample(obsB);'...
'obsPDist = ePDist?.tlim(tintObs);'...
'obsPitch = ePitch?.tlim(tintObs);'...
],ic)
zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsPDist = (obsPDist.time.epochUnix-mean(obsPDist.time.epochUnix))*CS_normal_velocity;
zObsPitch = (obsPitch.time.epochUnix-mean(obsPitch.time.epochUnix))*CS_normal_velocity;

% Model parameters
mms_2015Nov12.Bmodel;

%% Get pitchangles from mapped distributions

c_eval('ePitch?liou = tsFmap.pitchangles(mvaB?,15);',ic)
%%

h = irf_plot(4);
isub = 1;

elim = [10 1000];

if 1 % Magnetic field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');

  zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsB.data obsB.abs.data],'-');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,50)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Bx(0,0,zMod).^2+By(0,0,zMod).^2+Bz(0,0,zMod).^2)*1e9,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'B','(nT)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-13 13];
  end  
  if 1 % Electric field, no absolute value
    hca = h(isub); isub = isub + 1;
    % Colors
    B_colors = mms_colors('xyz1');

    set(hca,'colororder',B_colors)
    hca.ColorOrder = B_colors;
    linesObs = plot(hca,zObs,[obsE.data],'-');
    linesObs(1).Color = B_colors(1,:); 
    linesObs(2).Color = B_colors(2,:);
    linesObs(3).Color = B_colors(3,:);
    set(hca,'colororder',B_colors)
    %irf_legend(hca,{'E_L','E_M','E_N'},[0.01 0.2],'fontsize',14)
    irf_legend(hca,{'L','M','N'},[0.01 0.95],'fontsize',14)
    hold(hca,'on')
    zMod = linspace(-d*1.5,d*1.5,100)*3;
    %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
    lineMod = plot(hca,zMod*1e-3,[Ex(0,0,zMod)]*1e3,'--','color',B_colors(1,:));
    plot(hca,zMod*1e-3,[Ey(0,0,zMod)]*1e3,'--','color',B_colors(2,:))
    plot(hca,zMod*1e-3,[Ez(0,0,zMod)]*1e3,'--','color',B_colors(3,:))
    hold(hca,'off')
    %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
    %hca.Title.String = 'Magnetic field';
    hca.YGrid = 'on';
    hca.YLabel.String = {'E','(mV/m)'};
    hca.XLabel.String = 'N (km)';
    hca.XLim = [-30 30];
    hca.YLim = [-3.9 3.9];
  end  

  
if 1 % Distance: ePDist PSD pa low energies
  hca = h(isub); isub = isub + 1;    
  plotPitch = obsPitch.deflux.elim(elim).specrec('pa');
  pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))
  shading(hca,'flat')
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;

  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',12,'color',[0 0 0]);
  %hca.CLim = h(2).CLim;  
  xlabel(hca,'N (km)')  
end

if 1 % Distance: ePDist PSD pa low energies
  hca = h(isub); isub = isub + 1;    
  c_eval('plotPitch = ePitch?liou.deflux.elim(elim).specrec(''pa'');',ic)
  pcolor(hca,zf0,plotPitch.f,log10(plotPitch.p'))
  shading(hca,'flat')
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;

  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',12,'color',[0 0 0]);
  %hca.CLim = h(2).CLim;  
  xlabel(hca,'N (km)')  
end

h1pos = h(1).Position(3)*0.9;
for ip = 1:4    
  h(ip).Position(3) = h1pos*0.95;
  h(ip).XLim = [-100 100]; 
end

for ip = 3:4;
  %h(ip).YLim = [0 180];
  h(ip).XLim = [-100 100];
  h(ip).CLim = [7 8.1];
end

%% Single particle distributions
npanels = 4;
[h,h2] = initialize_combined_plot(npanels,2,2,0.4,'vertical');

isub = 1;

elim = [10 1000];

if 1 % Magnetic field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');

  zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsB.data obsB.abs.data],'-');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,50)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Bx(0,0,zMod).^2+By(0,0,zMod).^2+Bz(0,0,zMod).^2)*1e9,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'B','(nT)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-13 13];
  end  
if 1 % Electric field, no absolute value
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');

  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsE.data],'-');
  linesObs(1).Color = B_colors(1,:); 
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'E_L','E_M','E_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.01 0.95],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Ex(0,0,zMod)]*1e3,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[Ey(0,0,zMod)]*1e3,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Ez(0,0,zMod)]*1e3,'--','color',B_colors(3,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'E','(mV/m)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-3.9 3.9];
end  
if 1 % Distance: ePDist PSD pa low energies
  hca = h(isub); isub = isub + 1;    
  plotPitch = obsPitch.deflux.elim(elim).specrec('pa');
  pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))
  shading(hca,'flat')
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;

  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',12,'color',[0 0 0]);
  %hca.CLim = h(2).CLim;  
  xlabel(hca,'N (km)')  
end
if 1 % Distance: ePDist PSD pa low energies
  hca = h(isub); isub = isub + 1;    
  c_eval('plotPitch = ePitch?liou.deflux.elim(elim).specrec(''pa'');',ic)
  pcolor(hca,zf0,plotPitch.f,log10(plotPitch.p'))
  shading(hca,'flat')
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;

  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',12,'color',[0 0 0]);
  %hca.CLim = h(2).CLim;  
  xlabel(hca,'N (km)')  
end

h1pos = h(1).Position(3)*0.9;
for ip = 1:4    
  h(ip).Position(3) = h1pos*0.95;
  h(ip).XLim = [-100 100]; 
end
%%
for ip = 3:4;
  %h(ip).YLim = [0 180];
  h(ip).XLim = [-60 60];
  h(ip).CLim = [7 8.1];
end
%%
iPDS_toplot = 1;
isub = 1;
if 1 % Distance: ePDist PSD pa low energies
  hca = h2(isub); isub = isub + 1;    
  mms.plot_skymap(hca,tsFmap.deflux,'tint',times(1),'flat','energy',100,'vectors',{[1 0 0],'L';[0 1 0],'M';[0 0 1],'N'})
end
if 1 % Distance: ePDist PSD pa low energies
  hca = h2(isub); isub = isub + 1;    
  mms.plot_skymap(hca,tsFmap.deflux,'tint',times(1),'flat','energy',200,'vectors',{[1 0 0],'L';[0 1 0],'M';[0 0 1],'N'})
end
if 1 % Distance: ePDist PSD pa low energies
  hca = h2(isub); isub = isub + 1;    
  mms.plot_skymap(hca,obsPDist.deflux,'tint',times(1),'flat','energy',100,'vectors',{[1 0 0],'L';[0 1 0],'M';[0 0 1],'N'})
end
if 1 % Distance: ePDist PSD pa low energies
  hca = h2(isub); isub = isub + 1;    
  mms.plot_skymap(hca,obsPDist.deflux,'tint',times(1),'flat','energy',200,'vectors',{[1 0 0],'L';[0 1 0],'M';[0 0 1],'N'})
end