tint_boundary = irf.tint('2017-07-11T22:32:40.000Z/2017-07-11T22:33:50.000Z'); % LH separatrix
% Minvar of E
fhigh = 5;
flow = 0;
[out,l,v] = irf_minvar(gseE1.tlim(tint_boundary));
[outE,lE,vE] = irf_minvar(gseE1.tlim(tint_boundary).filt(flow,fhigh,[],5));
[outVi,lVi,vVi] = irf_minvar(gseE1.tlim(tint_boundary));

% I average B not B/|B|, to emphasize high B regions, dont know if this is
% optimal
Bav = mean(gseB1.tlim(tint_boundary).data,1); bav = Bav/norm(Bav);

%% Try out coordinate systems, make one common for entire boundary, easier for paper that way
newx = bav;
newy = cross(vE(1,:),bav);
newz = cross(newx,newy);
%newx = cross(newz,cross(v_direction,newz));
%newy = cross(newz,newx);
newxyz = [newx;newy;newz];
% vE does not give well separated min and intermediate direction, but decent normal direction
%newxyz = vVi;

legs = {'L','M','N'};

ic_tmp = ic;
ic = 1:4;
c_eval('bdryB? = gseB?*newxyz'';',ic);
c_eval('bdryE? = gseE?*newxyz'';',ic);
c_eval('bdryE?perp = gseE?perp*newxyz'';',ic);
c_eval('bdryE?par = gseE?par;',ic);
c_eval('bdryVExB? = gseVExB?*newxyz'';',ic);
c_eval('bdryVexB? = gseVexB?*newxyz'';',ic);
c_eval('bdryVe? = gseVe?*newxyz'';',ic);
c_eval('bdryVi? = gseVi?*newxyz'';',ic);
c_eval('bdryVe?perp = gseVe?perp*newxyz'';',ic);
c_eval('bdryVe?par = gseVe?par;',ic);
c_eval('bdryVi?perp = gseVi?perp*newxyz'';',ic);
c_eval('bdryE?perp = gseE?perp*newxyz'';',ic);
c_eval('bdryJ? = gseJ?*newxyz'';',ic);
bdryJcurl = gseJcurl*newxyz';
c_eval('bdryR? = gseR?*newxyz'';',ic);
c_eval('bdryRR? = gseRR?*newxyz'';',ic);

ic = ic_tmp;

ic = 1;
npanels = 6;
h = irf_plot(npanels);
tintZoom = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:38:00.00Z'); %20151112071854
tintZoom = tint_boundary;
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
  irf_legend(hca,legs,[0.98 0.9],'fontsize',12);  
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
if 0 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % beta
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
  irf_legend(hca,legs,[0.98 0.9],'fontsize',12);  
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
  irf_legend(hca,legs,[0.98 0.9],'fontsize',12);
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
if 1 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?perp.x.tlim(tint),bdryVe?perp.y.tlim(tint),bdryVe?perp.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,legs,[0.98 0.9],'fontsize',12);     
end
if 0 % Ve x B
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVExB?.x,bdryVExB?.y,bdryVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'-v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,legs,[0.98 0.9],'fontsize',12);     
end
if 1 % V ExB
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVExB?.x,bdryVExB?.y,bdryVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,legs,[0.98 0.9],'fontsize',12);     
end
if 1 % E perp
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryE?perp.x,bdryE?perp.y,bdryE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,legs,[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
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

