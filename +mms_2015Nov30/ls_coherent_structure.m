% E/B/phi lower hybrid analysis
tintZoom = irf.tint('2015-11-30T00:24:27.05Z/2015-11-30T00:24:27.20Z');
tintZoomThin = irf.tint('2015-11-30T00:24:27.10Z/2015-11-30T00:24:27.15Z');
tintZoomThin = irf.tint('2015-11-30T00:24:27.10Z/2015-11-30T00:24:27.15Z');
%tintZoomThin = irf.tint('2015-11-30T00:24:26.68Z/2015-11-30T00:24:26.75Z');
%tintZoomThin = irf.tint('2015-11-30T00:22:36.50Z/2015-11-30T00:22:37.50Z');
tintUTC = tintZoomThin.utc;

[phiEB,vbest,dirbest,thetas,corrs] = mms.lhwaveanalysis(tintZoomThin,gseE1,gseB1scm,gseB1,ne1,'plot',1,'lhfilt',10);

%%
tintZoom = tintZoomThin + 0.05*[-1 1];
x = dirbest;
z = mean(gseB1.tlim(tintZoomThin).data,1); z = z/norm(z);
y = cross(z,x);

vfactor = 1.5;
dcvPot1_ = irf.ts_scalar(dcvPot1.time,dcvPot1.data(:,1:4));

lmnE1 = gseE1*[x;y;z]'; lmnE1.name = 'fac E';

h=irf_plot({gseB1.tlim(tintZoom),...
            gseB1scm.tlim(tintZoom),...
            gseE1perp.tlim(tintZoom),...
            gseE1par.tlim(tintZoom),...
            lmnE1.tlim(tintZoom),...
            phiEB.*irf.ts_scalar(phiEB.time,repmat([vfactor 1],phiEB.length,1)),...
            -gseVexB1.tlim(tintZoom),...
            gseVe1.tlim(tintZoom),...
            gseVe1.tlim(tintZoom)+-dirbest*vbest*vfactor,...
            gseVe1.tlim(tintZoom)-gseVe1.resample(irf_time('2015-11-30T00:24:26.80Z','utc>epochtt')).data,...
            scPot1.tlim(tintZoom),...
            dcvPot1.tlim(tintZoom),...
            ne1.tlim(tintZoom)});
h(1).YLabel.String ='B';          
h(2).YLabel.String ='scm B';          
irf_legend(h(5),{'Ek','E norm','Epar'},[0.95 0.05])
h(6).YLabel.String ='Phi'; irf_legend(h(6),{'phi E','phi B'},[0.95 0.05])
h(7).YLabel.String ='-vexB';
h(8).YLabel.String ='ve';
h(9).YLabel.String ='ve-vph';
h(10).YLabel.String ='ve-veref';
h(11).YLabel.String ='scPot';

add_length_on_top(h(1),vbest*vfactor,1)

%% Figure
% 2-3 panels overview, slightly longer time interval
tintZoom = irf.tint('2015-11-30T00:24:20.00Z/2015-11-30T00:24:30.00Z');
tintZoomThin = irf.tint('2015-11-30T00:24:27.10Z/2015-11-30T00:24:27.15Z');
tintZoomPhi = irf.tint('2015-11-30T00:24:20.00Z/2015-11-30T00:24:30.00Z');
vfactor = 1.5;

x = dirbest;
z = mean(gseB1.tlim(tintZoomThin).data,1); z = z/norm(z);
y = cross(z,x);

lmnE1 = gseE1*[x;y;z]'; lmnE1.name = 'fac E';
lmnB1 = gseB1scm*[x;y;z]'; lmnB1.name = 'fac B';


npanels = 9;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
if 1 % B
  iisub = iisub+1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?.x.tlim(tintZoom),gseB?.y.tlim(tintZoom),gseB?.z.tlim(tintZoom)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  iisub = iisub+1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('xy'))
  c_eval('irf_plot(hca,{ne?.tlim(tintZoom)},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('xy'))    
end
if 0 % Vi
  iisub = iisub+1;
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve
  iisub = iisub+1;
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)
  hca.YLabel.String = {'V_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % J ExB
  iisub = iisub+1;
  c_eval('gseJExB?filt = gseJExB?.filt(0,100,[],5)',ic)
  hca = irf_panel('J ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJExB?filt.x.tlim(tintZoom)*1e-3,gseJExB?filt.y.tlim(tintZoom)*1e-3,gseJExB?filt.z.tlim(tintZoom)*1e-3},''comp'');',ic)
  %c_eval('hExB = irf_plot(hca,{gseVExB?.x.tlim(tintZoom),gseVExB?.y.tlim(tintZoom),gseVExB?.z.tlim(tintZoom)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{dbcsVe?.x,dbcsVe?.y,dbcsVe?.z},''comp'');',ic)
  hca.YLabel.String = {'J_{e,ExB}','(\mu A/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hh = irf_legend(hca,{['f  <100 Hz']},[0.98 0.1],'fontsize',12); hh.Color = [0 0 0]; 
end

hca = irf_panel('empty');
ffilt = 3;

if 1 % B scm  
  c_eval('facBfilt = lmnB?.filt(ffilt,0,[],5);',ic)
  hca = irf_panel('scm B');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  %c_eval('irf_plot(hca,{gseB?scm.tlim(tintZoomThin).x,gseB?scm.tlim(tintZoomThin).y,gseB?scm.tlim(tintZoomThin).z},''comp'');',ic)
  c_eval('irf_plot(hca,{facBfilt.tlim(tintZoomThin).x,facBfilt.tlim(tintZoomThin).y,facBfilt.tlim(tintZoomThin).z},''comp'');',ic)
  hca.YLabel.String = {'B_{SCM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'k','||\times k','||'},[0.98 0.9],'fontsize',12);
  hh = irf_legend(hca,{['f  > ' num2str(ffilt) ' Hz']},[0.1 0.9],'fontsize',12); hh.Color = [0 0 0];
end
if 1 % E
  c_eval('facEfilt = lmnE?.filt(ffilt,0,[],5);',ic)
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseE?.tlim(tintZoomThin).x,gseE?.tlim(tintZoomThin).y,gseE?.tlim(tintZoomThin).z},''comp'');',ic)
  c_eval('irf_plot(hca,{facEfilt.tlim(tintZoomThin).x,facEfilt.tlim(tintZoomThin).y,facEfilt.tlim(tintZoomThin).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'\perp,1','\perp,2','||'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'k','||\times k','||'},[0.98 0.6],'fontsize',12);
  hh = irf_legend(hca,{['f  > ' num2str(ffilt) ' Hz']},[0.1 0.9],'fontsize',12); hh.Color = [0 0 0];
end
if 1 % phi
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('12'))
  %c_eval('irf_plot(hca,{gseE?.tlim(tintZoomThin).x,gseE?.tlim(tintZoomThin).y,gseE?.tlim(tintZoomThin).z},''comp'');',ic)
  phiE = phiEB.clone(phiEB.time,phiEB.data(:,1)*vfactor);
  phiB = phiEB.clone(phiEB.time,phiEB.data(:,2));
  irf_plot(hca,{phiE,phiB},'comp');
  ylabel(hca,{'\phi','(V)'},'interpreter','tex')
  %hca.YLabel.String = {'\phi','(V)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'\phi_E','\phi_B'},[0.98 0.5],'fontsize',12);
  irf_legend(hca,{['v = ' num2str(vbest*vfactor,'%.0f') ' km/s \times [' num2str(dirbest(1),'%.1f') ' '  num2str(dirbest(2),'%.1f') ' '  num2str(dirbest(3),'%.1f') ']']},[0.95 0.95],'fontsize',12);
end
if 1 % sc Pot
  hca = irf_panel('scPot');
  set(hca,'ColorOrder',mms_colors('xy'))
  c_eval('irf_plot(hca,{scPot?},''comp'');',ic)
  hca.YLabel.String = {'sc Pot','(V)'};
  set(hca,'ColorOrder',mms_colors('xy'))    
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:iisub iisub+2:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

for  ii = 1:npanels
  h(ii).XGrid = 'off';
  h(ii).YGrid = 'off';
end

irf_zoom(h(1:iisub),'x',tintZoom)
irf_zoom(h,'y')

irf_pl_mark(h(1:iisub),tintZoomThin,'g')

irf_zoom(h(iisub+2:npanels),'x',tintZoomThin)

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))
hca = irf_panel('empty');  hca.Visible = 'off';

h_length = add_length_on_top(h(iisub+2),vbest*vfactor,0.5);
h_length.XTickLabel{end} = ['     ' h_length.XTickLabel{end} ' km'];
irf_zoom(h,'y')

%% Estimate half widths
tintZoomThinner = irf.tint('2015-11-30T00:24:27.116Z/2015-11-30T00:24:27.14Z');
phiE = phiEB.clone(phiEB.time,phiEB.data(:,1)*vfactor).tlim(tintZoomThinner);
phiB = phiEB.clone(phiEB.time,phiEB.data(:,2)).tlim(tintZoomThinner);

t0 = irf_time('2015-11-30T00:24:27.126Z','utc>epochtt');
r = (phiE.time-t0)*vbest*vfactor;
z = 0;

phi = @(phi0,r,Lr,z,Lz) phi0*exp(-r.^2/2/Lr^2-z.^2/2/Lz^2);
Lr = 2; % km
Lz = 20; % km
phi0 = -15;
phiFit = phi(phi0,r,Lr,z,Lz);
tsPhiFit = irf.ts_scalar(phiE.time,phiFit);

% density;
n = @(phi0,r,Lr,z,Lz) units.eps0*phi0/Lr/Lr*(1-r.^2/Lr/Lr)*exp(-r.^2/2/Lr^2-z.^2/2/Lz^2);
nFit = phi(phi0,r,Lr,z,Lz);
tsNFit = irf.ts_scalar(phiE.time,nFit);


h = irf_plot(2);
hca = irf_panel('phi');
irf_plot(hca,{phiE,phiB,tsPhiFit},'comp');

hca = irf_panel('n');
irf_plot(hca,{tsNFit},'comp');

