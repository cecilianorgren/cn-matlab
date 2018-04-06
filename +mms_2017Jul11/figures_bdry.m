%% Comparison between ion and electron pressure
nrows = 2;
ncols = 3;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end

timeline = ni1.time;

isub = 1;

if 1 % PB vs Pi
  hca = h(isub); isub = isub + 1;
  scatter(hca,PB1.resample(timeline).data,gsePi1.trace.resample(timeline).data/3,'.')
  hca.XLabel.String = 'P_B';
  hca.YLabel.String = 'P_i';
  hold(hca,'on')
  beta0 = 0.16;
  betaslope = 1;
  hl = plot(hca,hca.XLim,beta0-betaslope*hca.XLim,'linewidth',2);
  hl2 = plot(hca,0.011*[1 1],hca.YLim,'linewidth',2);
  hold(hca,'off')
end
if 1 % PB vs Pe
  hca = h(isub); isub = isub + 1;
  scatter(hca,PB1.resample(timeline).data,gsePe1.trace.resample(timeline).data/3,'.')
  hca.XLabel.String = 'P_B';
  hca.YLabel.String = 'P_2';
  hold(hca,'on')
  beta0 = 0.04;
  betaslope = 0.3;
  hl = plot(hca,hca.XLim,beta0-betaslope*hca.XLim,'linewidth',2);
  hl2 = plot(hca,0.011*[1 1],hca.YLim,'linewidth',2);
  hold(hca,'off')
end
if 1 % Pe vs Pi with colorcoding according to PB
  hca = h(isub); isub = isub + 1;
  PBlim = 0.012;
  plotx = gsePe1.trace.resample(timeline).data/3;
  ploty = gsePi1.trace.resample(timeline).data/3;
  plots = plotx*0+1;
  plotc = PB1.resample(timeline).data;
  if 1
    plotc(plotc<PBlim) = 0;
    plotc(plotc>PBlim) = 0.8;
    clim = [0 1];
  end
  %plotc(plotc<PBlim)= NaN;
  scatter(hca,plotx,ploty,plots,plotc,'.')
  hca.CLim = clim;
  hca.XLabel.String = 'P_e';
  hca.YLabel.String = 'P_i';
  hcb =colorbar('peer',hca);
  hcb.Position = hcb.Position;
end
if 1
  hca = h(isub); isub = isub + 1;
  scatter(hca,gsePe1.trace.resample(timeline).data/3,gsePi1.trace.resample(timeline).data/3,'.')
  hca.XLabel.String = 'P_e';
  hca.YLabel.String = 'P_i';
  %hca.XLim = hca.YLim;
end
if 0
  hca = h(isub); isub = isub + 1;
  scatter(hca,gseTe1.trace.resample(timeline).data/3,gseTi1.resample(timeline).trace.data/3)
  hca.XLabel.String = 'T_e';
  hca.YLabel.String = 'T_i';
  hca.XLim = hca.YLim;
end
if 0
  hca = h(isub); isub = isub + 1;
  scatter(hca,ne1.resample(timeline).data/3,ni1.resample(timeline).data/3)
  hca.XLabel.String = 'n_e';
  hca.YLabel.String = 'n_i';
  hca.XLim = hca.YLim;
end
if 0
  hca = h(isub); isub = isub + 1;
  scatter(hca,ni1.resample(timeline).data/3,gsePi1.resample(timeline).trace.data/3)
  hca.XLabel.String = 'n_i';
  hca.YLabel.String = 'P_i';
end
if 0
  hca = h(isub); isub = isub + 1;
  scatter3(hca,ni1.resample(timeline).data/3,gsePi1.resampleh(timeline).trace.data/3,gseTi1.resample(timeline).trace.data/3)
  hca.XLabel.String = 'n_i';
  hca.YLabel.String = 'P_i';
  hca.ZLabel.String = 'T_i';
end
for ip = 1:npanels
  axis(h(ip),'square')
end

%% Pressure balance
ic = 1;
npanels = 4;
h = irf_plot(npanels);
tintZoom = irf.tint('2017-07-11T22:33:20.00Z/2017-07-11T22:33:45.00Z'); %20151112071854
tint = tintZoom;

if 1 % iPDist deflux omni
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  c_eval('irf_plot(hca,gseTi?.trace/3,''k'');',ic)
  hold(hca,'off')
end
if 1 % ePDist deflux omni
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet') 
  hold(hca,'on')
  c_eval('irf_plot(hca,gseTe?.trace/3,''k'');',ic)
  hold(hca,'off')
end
if 1 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{beta1.tlim(tint),beta1i.tlim(tint),beta1e.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('\beta',ic),'(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'i+e','i','e'},[0.98 0.9],'fontsize',12);
end
if 1 % Pressures
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{gsePe?.trace/3,facPe?.xx,0.5*(facPe?.yy+facPe?.zz),'...
                         'gsePi?.trace/3,PB?,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_e','P_{e||}','P_{e\perp}','P_i','P_B','P_{tot}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
end
if 0 % Pressures, stacked
  hca = irf_panel('Pressure stacked');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_patch(hca,{gsePe?.trace/3,facPe?.xx,0.5*(facPe?.yy+facPe?.zz),'...
                         'gsePi?.trace/3,PB?,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3});'],ic)
  irf_patch(hca,{PB1,gsePi1.trace/3,gsePe1.trace/3})
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_e','P_{e||}','P_{e\perp}','P_i','P_B','P_{tot}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
end
if 0 % Pressures, total, 4sc
  hca = irf_panel('Pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePe1.trace/3+PB1.resample(gsePe1)+gsePi1.resample(gsePe1).trace/3,...    
                gsePe2.trace/3+PB2.resample(gsePe2)+gsePi1.resample(gsePe2).trace/3,...    
                gsePe3.trace/3+PB3.resample(gsePe3)+gsePi1.resample(gsePe3).trace/3,...    
                gsePe4.trace/3+PB4.resample(gsePe4)+gsePi1.resample(gsePe4).trace/3},'comp');
  hca.YLabel.String = {'P_e+P_i+P_B','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Ion pressures, total, 4sc
  hca = irf_panel('Ion pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePi1.trace/3,gsePi2.trace/3,gsePi3.trace/3,gsePi4.trace/3},'comp');
  hca.YLabel.String = {'P_i','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Electron pressures, total, 4sc
  hca = irf_panel('Pe, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePe1.trace/3,gsePe2.trace/3,gsePe3.trace/3,gsePe4.trace/3},'comp');
  hca.YLabel.String = {'P_e','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Electron temperature, total, 4sc
  hca = irf_panel('Te, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseTe1.trace/3,gseTe2.trace/3,gseTe3.trace/3,gseTe4.trace/3},'comp');
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % gradPi
  hca = irf_panel('gradPi');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPi.x*1e3,gseGradPi.y*1e3,gseGradPi.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_i','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % E
  hca = irf_panel('E lowf');
  ffiltE = 0.5;
  c_eval('gseE?highf = gseE?.filt(3,0,[],3);',ic)
  c_eval('gseE?lowf = gseE?-gseE?highf;',ic)
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?lowf.x,gseE?lowf.y,gseE?lowf.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f<%g Hz',ffiltE),[0.98 0.1],'fontsize',12);
end

irf_zoom(h,'x',tintZoom)
%irf_zoom(h,'y')
irf_plot_axis_align

h(1).Title.String = {['MMS' num2str(ic)],''};

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
for ii = 1:npanels
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0])
  h(ii).FontSize = 14;  
  h(ii).YLabel.FontSize = 13;
end

%% Gradient from magnetic field
npanels = 7;
h = irf_plot(npanels);

if 1 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.abs.tlim(tint),gseB2.abs.tlim(tint),gseB3.abs.tlim(tint),gseB4.abs.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{z}','(nT)'};
end
if 1 % BX
  hca = irf_panel('BX');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.x.tlim(tint),gseB2.x.tlim(tint),gseB3.x.tlim(tint),gseB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{x}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BY
  hca = irf_panel('BY');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.y.tlim(tint),gseB2.y.tlim(tint),gseB3.y.tlim(tint),gseB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{y}','(nT)'};
end
if 1 % BZ
  hca = irf_panel('BZ');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.z.tlim(tint),gseB2.z.tlim(tint),gseB3.z.tlim(tint),gseB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{z}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % JM
  hca = irf_panel('JM');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{gseJ1.y,gseJ2.y,gseJ3.y,gseJ4.y,gseAvJ.y,gseJcurl.y},'comp');   
  lines = irf_plot(hca,{gseJ1.y,gseJ2.y,gseJ3.y,gseJ4.y},'comp');   
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1par.tlim(tint),gseVe2par.tlim(tint),gseVe3par.tlim(tint),gseVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp X
  hca = irf_panel('Ve perp X');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1perp.tlim(tint).x,gseVe2perp.tlim(tint).x,gseVe3perp.tlim(tint).x,gseVe4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,x}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp Y
  hca = irf_panel('Ve perp Y');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1perp.tlim(tint).y,gseVe2perp.tlim(tint).y,gseVe3perp.tlim(tint).y,gseVe4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,y}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp Z
  hca = irf_panel('Ve perp Z');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1perp.tlim(tint).z,gseVe2perp.tlim(tint).z,gseVe3perp.tlim(tint).z,gseVe4perp.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,z}','(km/s)'},'interpreter','tex');
end
if 0 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
 irf_plot(hca,{gseEVexB1.abs,gseEVexB2.abs,gseEVexB3.abs,gseEVexB4.abs},'comp'); 
  hca.YLabel.String = {'|E+v_{e}xB|','(mV/m)'};
  %hca.YLabel.String = {'E''_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234')) 
  %hca.YLim = [-10 10];
end
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end
if 1 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  lines = irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');   
  hca.YLabel.String = {'J_{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P',['(10^{-3} ' gseGradPe.units ')']};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % scPot
  hca = irf_panel('scPot');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-1*scPot1.tlim(tint),-1*scPot2.tlim(tint),-1*scPot3.tlim(tint),-1*scPot4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('-scPot',ic),'(V)'};
end


irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Plasma Beta
npanels = 7;
h = irf_plot(npanels);

if 1 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.abs.tlim(tint),gseB2.abs.tlim(tint),gseB3.abs.tlim(tint),gseB4.abs.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{z}','(nT)'};
end
if 1 % BX
  hca = irf_panel('BX');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.x.tlim(tint),gseB2.x.tlim(tint),gseB3.x.tlim(tint),gseB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{x}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BY
  hca = irf_panel('BY');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.y.tlim(tint),gseB2.y.tlim(tint),gseB3.y.tlim(tint),gseB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{y}','(nT)'};
end
if 1 % BZ
  hca = irf_panel('BZ');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.z.tlim(tint),gseB2.z.tlim(tint),gseB3.z.tlim(tint),gseB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{z}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % JM
  hca = irf_panel('JM');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{gseJ1.y,gseJ2.y,gseJ3.y,gseJ4.y,gseAvJ.y,gseJcurl.y},'comp');   
  lines = irf_plot(hca,{gseJ1.y,gseJ2.y,gseJ3.y,gseJ4.y},'comp');   
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1par.tlim(tint),gseVe2par.tlim(tint),gseVe3par.tlim(tint),gseVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp X
  hca = irf_panel('Ve perp X');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1perp.tlim(tint).x,gseVe2perp.tlim(tint).x,gseVe3perp.tlim(tint).x,gseVe4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,x}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp Y
  hca = irf_panel('Ve perp Y');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1perp.tlim(tint).y,gseVe2perp.tlim(tint).y,gseVe3perp.tlim(tint).y,gseVe4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,y}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp Z
  hca = irf_panel('Ve perp Z');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1perp.tlim(tint).z,gseVe2perp.tlim(tint).z,gseVe3perp.tlim(tint).z,gseVe4perp.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,z}','(km/s)'},'interpreter','tex');
end
if 0 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
 irf_plot(hca,{gseEVexB1.abs,gseEVexB2.abs,gseEVexB3.abs,gseEVexB4.abs},'comp'); 
  hca.YLabel.String = {'|E+v_{e}xB|','(mV/m)'};
  %hca.YLabel.String = {'E''_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234')) 
  %hca.YLim = [-10 10];
end
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end
if 1 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  lines = irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');   
  hca.YLabel.String = {'J_{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P',['(10^{-3} ' gseGradPe.units ')']};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % scPot
  hca = irf_panel('scPot');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-1*scPot1.tlim(tint),-1*scPot2.tlim(tint),-1*scPot3.tlim(tint),-1*scPot4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('-scPot',ic),'(V)'};
end


irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Spacecraft position in local coordinate system
% Subplot spacecraft position
nrows = 1;
ncols = 3;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end

isub = 1;

xyzlim = ceil((max([abs(bdryRR1), abs(bdryRR2), abs(bdryRR3), abs(bdryRR4)]))/5)*5;
markersize = 12;
linewidth = 1.5;
fontsize = 15;

hca = h(isub); isub = isub + 1;
hold(hca,'on')
plot(hca,mean(bdryRR1(1)),mean(bdryRR1(2)),'ks','markersize',markersize,'linewidth',linewidth)
plot(hca,mean(bdryRR2(1)),mean(bdryRR2(2)),'rd','markersize',markersize,'linewidth',linewidth)
plot(hca,mean(bdryRR3(1)),mean(bdryRR3(2)),'go','markersize',markersize,'linewidth',linewidth)
plot(hca,mean(bdryRR4(1)),mean(bdryRR4(2)),'bv','markersize',markersize,'linewidth',linewidth)
xlabel(hca,'x / v')
ylabel(hca,'y / N')

hca = h(isub); isub = isub + 1;
hold(hca,'on')
plot(hca,mean(bdryRR1(1)),mean(bdryRR1(3)),'ks','markersize',markersize,'linewidth',linewidth)
plot(hca,mean(bdryRR2(1)),mean(bdryRR2(3)),'rd','markersize',markersize,'linewidth',linewidth)
plot(hca,mean(bdryRR3(1)),mean(bdryRR3(3)),'go','markersize',markersize,'linewidth',linewidth)
plot(hca,mean(bdryRR4(1)),mean(bdryRR4(3)),'bv','markersize',markersize,'linewidth',linewidth)
xlabel(hca,'x / v')
ylabel(hca,'z / B')

hca = h(isub); isub = isub + 1;
hold(hca,'on')
plot(hca,mean(bdryRR1(2)),mean(bdryRR1(3)),'ks','markersize',markersize,'linewidth',linewidth)
plot(hca,mean(bdryRR2(2)),mean(bdryRR2(3)),'rd','markersize',markersize,'linewidth',linewidth)
plot(hca,mean(bdryRR3(2)),mean(bdryRR3(3)),'go','markersize',markersize,'linewidth',linewidth)
plot(hca,mean(bdryRR4(2)),mean(bdryRR4(3)),'bv','markersize',markersize,'linewidth',linewidth)
xlabel(hca,'y / N')
ylabel(hca,'z / B')

for ip = 1:npanels
  hca = h(ip);
  box(hca,'on')
  axis(hca,'square')
  hca.XLim = xyzlim*[-1 1];
  hca.YLim = xyzlim*[-1 1];
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.FontSize = fontsize;
end

%% Plasma parameters, beta, pressure etc
ic = 1;
npanels = 4;
h = irf_plot(npanels);
tintZoom = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:38:00.00Z'); %20151112071854
tint = tintZoom;

zoomy = [];
isub = 0;
if 1 % B  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryB?.x.tlim(tint),bdryB?.y.tlim(tint),bdryB?.z.tlim(tint),bdryB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B','|B|'},[0.98 0.9],'fontsize',12);  
end
if 0 % iPDist deflux omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  c_eval('irf_plot(hca,gseTi?.trace/3,''k'');',ic)
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};  
end
if 0 % ePDist deflux omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet') 
  hold(hca,'on')
  c_eval('irf_plot(hca,gseTe?.trace/3,''k'');',ic)
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};  
end
if 0 % ne ni
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))  
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
end
if 1 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % beta
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{beta1.tlim(tint),beta1i.tlim(tint),beta1e.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('\beta',ic),'(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'i+e','i','e'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log';
  hca.YLim = [1e-2 1e3];
  hca.YTick = 10.^[-1:2];
end
if 0 % Pressures
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 0 % Pressures, stacked
  hca = irf_panel('Pressure stacked');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_patch(hca,{gsePe?.trace/3,facPe?.xx,0.5*(facPe?.yy+facPe?.zz),'...
                         'gsePi?.trace/3,PB?,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3});'],ic)
  irf_patch(hca,{PB1,gsePi1.trace/3,gsePe1.trace/3})
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_e','P_{e||}','P_{e\perp}','P_i','P_B','P_{tot}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
end
if 0 % Pressures, total, 4sc
  hca = irf_panel('Pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePe1.trace/3+PB1.resample(gsePe1)+gsePi1.resample(gsePe1).trace/3,...    
                gsePe2.trace/3+PB2.resample(gsePe2)+gsePi1.resample(gsePe2).trace/3,...    
                gsePe3.trace/3+PB3.resample(gsePe3)+gsePi1.resample(gsePe3).trace/3,...    
                gsePe4.trace/3+PB4.resample(gsePe4)+gsePi1.resample(gsePe4).trace/3},'comp');
  hca.YLabel.String = {'P_e+P_i+P_B','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Ion pressures, total, 4sc
  hca = irf_panel('Ion pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePi1.trace/3,gsePi2.trace/3,gsePi3.trace/3,gsePi4.trace/3},'comp');
  hca.YLabel.String = {'P_i','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Electron pressures, total, 4sc
  hca = irf_panel('Pe, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePe1.trace/3,gsePe2.trace/3,gsePe3.trace/3,gsePe4.trace/3},'comp');
  hca.YLabel.String = {'P_e','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Electron temperature, total, 4sc
  hca = irf_panel('Te, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseTe1.trace/3,gseTe2.trace/3,gseTe3.trace/3,gseTe4.trace/3},'comp');
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % gradPi
  hca = irf_panel('gradPi');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPi.x*1e3,gseGradPi.y*1e3,gseGradPi.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_i','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryE?.x,bdryE?.y,bdryE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);  
end
if 0 % E
  hca = irf_panel('E lowf');
  ffiltE = 0.5;
  c_eval('gseE?highf = gseE?.filt(3,0,[],3);',ic)
  c_eval('gseE?lowf = gseE?-gseE?highf;',ic)
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?lowf.x,gseE?lowf.y,gseE?lowf.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f<%g Hz',ffiltE),[0.98 0.1],'fontsize',12);
end
if 0 % Lp from Eperp
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('Lp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Lp?*1e-3},''comp'');',ic)
  hca.YLabel.String = {'L_p','(km)'};  
  hca.YScale = 'log';
  hca.YLim = [1e2 5e3];
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
irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align

h(1).Title.String = {['MMS' num2str(ic)],''};

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
for ii = 1:npanels
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0])
  h(ii).FontSize = 14;  
  h(ii).YLabel.FontSize = 13;
end

%% Velocities in new coordinate system, bdry
npanels = 7;
h = irf_plot(npanels);
ic = 1;

tint = irf.tint('2017-07-11T22:33:15.00Z/2017-07-11T22:33:35.00Z'); %20151112071854

pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{bdryB?.x,bdryB?.y,bdryB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);
end 
if 1 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{bdryJcurl.x,bdryJcurl.y,bdryJcurl.z},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 1 % J  
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryJ?.x,bdryJ?.y,bdryJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi  
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVi?.x,bdryVi?.y,bdryVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12); 
end
if 1 % Vi  
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVi?perp.x,bdryVi?perp.y,bdryVi?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12); 
end
if 0 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?.x.tlim(tint),bdryVe?.y.tlim(tint),bdryVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);     
end
if 1 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?perp.x.tlim(tint),bdryVe?perp.y.tlim(tint),bdryVe?perp.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve x B
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVExB?.x,bdryVExB?.y,bdryVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'-v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);     
end
if 1 % V ExB
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVExB?.x,bdryVExB?.y,bdryVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);     
end
if 1 % E perp
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryE?perp.x,bdryE?perp.y,bdryE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','||'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end


irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Velocities in new coordinate system, mva
npanels = 8;
h = irf_plot(npanels);
ic = 1;

tint = irf.tint('2017-07-11T22:33:15.00Z/2017-07-11T22:33:35.00Z'); %20151112071854
tint = irf.tint('2017-07-11T22:33:15.00Z/2017-07-11T22:33:35.00Z'); %20151112071854

pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 1 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{mvaB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('1'))
end 
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end 
if 1 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{mvaJcurl.x,mvaJcurl.y,mvaJcurl.z},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 1 % J  
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaJ?.x,mvaJ?.y,mvaJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi  
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x,mvaVi?.y,mvaVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi  
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?perp.x,mvaVi?perp.y,mvaVi?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x,mvaVe?perp.y,mvaVe?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve x B
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x,mvaVExB?.y,mvaVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'-v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % V ExB
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x,mvaVExB?.y,mvaVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);   
end
if 1 % E perp
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?perp.x,mvaE?perp.y,mvaE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end


irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Plot overview figure with focus on electrons, including single time electron distributions
npanels = 10;
tint = irf.tint('2017-07-11T22:33:00.00Z/2017-07-11T22:33:35.00Z'); %20151112071854
tintZoom = irf.tint('2017-07-11T22:31:50.00Z/2017-07-11T22:33:35.20Z'); % LH separatrix

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
  c_eval('irf_plot(hca,{bdryB?.x.tlim(tint),bdryB?.y.tlim(tint),bdryB?.z.tlim(tint),bdryB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B','|B|'},[0.98 0.9],'fontsize',12);  
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
if 1 % E perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryE?perp.x,bdryE?perp.y,bdryE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 1 % J  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryJ?.x.tlim(tint),bdryJ?.y.tlim(tint),bdryJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);    
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
if 1 % ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?perp.x.tlim(tint),bdryVe?perp.y.tlim(tint),bdryVe?perp.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);   
end
if 1 % Ve par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?par.tlim(tint)},''comp'');',ic)
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
if 1 % ePDist pa 64
  isub = isub + 1;
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [100 30000];  
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % Te par perp Ti/Tref
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
tintDist = irf.tint('2017-07-11T22:33:19.00Z/2017-07-11T22:33:23.00Z');
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('dist_scm = ePDist?;',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
c_eval('ePitch = ePitch?.convertto(''s^3/km^6'');',ic)

% Plot format input
vlim = 50*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [-3 1.5];  
projclim_int = [-13 -9];
  

c_eval('times = ePDist?.time;',ic)
[tind,~] = times.tlim(tintDist);

doReducedF = 1;
for it_ = 1%:numel(tind)
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
  perp1 = cross(par,newy);
  perp2 = cross(par,perp1);  
  

  timeUTC = time.utc;      
  isub = 1;

  %vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B'};
  vectors = {[1 0 0],'x'; [0 1 0],'y'; [0 0 1],'z'};
  
  if 1 % Perpendicular plane, slice
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{BxN}','v_{Bx(BxN)}','v_{||}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);
    hca.Title.String = '';
  end 
  if 1 % B plane 1, slice
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{BxN}','v_{||}','-v_{Bx(BxN)}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
    hca.Title.String = '';
  end 
  if 1 % B plane 2, slice
    hca = h2(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{Bx(BxN)}','v_{||}','v_{BxN}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
    hca.Title.String = '';
  end
  if 1 % Perpendicular plane, integrated
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{BxN}','v_{Bx(BxN)}','v_{||}'};
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = '';
  end 
  if 1 % B plane 1, integrated
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{BxN}','v_{||}','-v_{Bx(BxN)}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = '';
  end  
  if 1 % B plane 2, integrated
    hca = h2(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{Bx(BxN)}','v_{||}','v_{BxN}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = '';
  end  
  if 1 % Perpendicular plane, integrated, smaller vint
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{BxN}','v_{Bx(BxN)}','v_{||}'};
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = '';
  end 
  vint = 10000*[-1 1];
  if 1 % B plane 1, integrated, smaller vint
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{BxN}','v_{||}','-v_{Bx(BxN)}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);
    hca.Title.String = '';
  end  
  if 1 % B plane 2, integrated, smaller vint
    hca = h2(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{Bx(BxN)}','v_{||}','v_{BxN}'};
    %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    mms.plot_int_projection(hca,dist_scm,'t',time,'xyz',xyz,'vlim',vlim,'vzint',vint,'clim',projclim_int,'vlabel',vlabels,'colorbar',1);        
    hca.Title.String = '';
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
  
  h2(1).Title.String = timeUTC(1:23);
  for ii = 1:4
    colormap(h2(ii),strCMap)
  end
  cn.print(['e_proj_in_fix_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','path',[eventPath 'proj_int_mms1/'])
end

%% Earliest package of LH waves
npanels = 7;
h = irf_plot(npanels);
ic = 1;

tint = irf.tint('2017-07-11T22:33:00.00Z/2017-07-11T22:33:25.00Z'); %20151112071854

pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{bdryB?.x,bdryB?.y,bdryB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);
end 
if 0 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{bdryJcurl.x,bdryJcurl.y,bdryJcurl.z},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 1 % J  
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryJ?.x,bdryJ?.y,bdryJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi  
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVi?.x,bdryVi?.y,bdryVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12); 
end
if 0 % Vi  
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVi?perp.x,bdryVi?perp.y,bdryVi?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12); 
end
if 0 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?.x.tlim(tint),bdryVe?.y.tlim(tint),bdryVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?perp.x.tlim(tint),bdryVe?perp.y.tlim(tint),bdryVe?perp.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve x B
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVExB?.x,bdryVExB?.y,bdryVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'-v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);     
end
if 0 % V ExB
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVExB?.x,bdryVExB?.y,bdryVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);     
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % E perp
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryE?perp.x,bdryE?perp.y,bdryE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','||'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E par, low freq
  hca = irf_panel('E par, low freq');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 10;
  c_eval('irf_plot(hca,{bdryE?par.filt(0,ffilt,[],3)},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{sprintf('f<%g Hz',ffilt)},[0.05 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E par, high freq
  hca = irf_panel('E par, high freq');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 10;
  c_eval('irf_plot(hca,{bdryE?par.filt(ffilt,0,[],3)},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{sprintf('f>%g Hz',ffilt)},[0.05 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(km/s)'};
end


irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end
