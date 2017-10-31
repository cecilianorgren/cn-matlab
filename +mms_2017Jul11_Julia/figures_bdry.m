%% Pressure balance
ic = 1;
npanels = 6;
h = irf_plot(npanels);
tintZoom = irf.tint('2017-07-11T22:33:28.00Z/2017-07-11T22:34:10.00Z'); %20151112071854


if 1 % iPDist deflux omni
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet') 
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
if 1 % Ion pressures, total, 4sc
  hca = irf_panel('Ion pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePi1.trace/3,gsePi2.trace/3,gsePi3.trace/3,gsePi4.trace/3},'comp');
  hca.YLabel.String = {'P_i','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % Electron pressures, total, 4sc
  hca = irf_panel('Pe, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePe1.trace/3,gsePe2.trace/3,gsePe3.trace/3,gsePe4.trace/3},'comp');
  hca.YLabel.String = {'P_e','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % Electron temperature, total, 4sc
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

%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

