%% US-Japan

% Rotate things into new coordinate system 
L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9239];
lmn = [L;M;N];


% Rotate data
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)

c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';',ic)
c_eval('mvaVe? = gseVe?*lmn''; mvaVe?.name = ''Ve LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
c_eval('mvaVExB? = gseVExB?*lmn''; mvaVExB?.name = ''VExB LMN'';',ic)
%c_eval('mvaJ? = gseJ?*lmn''; mvaJ?.name = ''J LMN'';',ic)
%c_eval('mvaJe? = gseJe?*lmn''; mvaJe?.name = ''Je LMN'';',ic)
%c_eval('mvaJi? = gseJi?*lmn''; mvaJi?.name = ''Ji LMN'';',ic)
%mvaJcurl = gseJcurl*lmn'; mvaJcurl.name = 'J LMN CURL';
c_eval('mvaPi? = lmn*gsePi?*lmn''; mvaPi?.units = gsePi?.units;',ic)
c_eval('mvaPe? = lmn*gsePe?*lmn''; mvaPe?.units = gsePe?.units;',ic)
c_eval('mvaTi? = lmn*gseTi?*lmn''; mvaTi?.units = gseTi?.units;',ic)
c_eval('mvaTe? = lmn*gseTe?*lmn''; mvaTe?.units = gseTe?.units;',ic)

c_eval('mvaPe? = lmn*gsePe?*lmn''; mvaPe?.units = gsePe?.units;',ic)


c_eval('mvaR? = gseR?*lmn''; mvaR?.name = ''R LMN'';',ic)

%% Figure: Overview 1
ic = 1;

tint_edr = irf.tint('2017-07-11T22:33:58.00Z/2017-07-11T22:34:10.00Z'); %20151112071854

npanels = 4;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
fontsize = 16;

isub = 0;
zoomy = [];

if 1 % B LMN
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize);
end 
if 1 % Ve
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize);
end
if 1 % Pe-off LMN
  hca = irf_panel('Pe off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPe?.xy.tlim(tint)*1e3,mvaPe?.xz.tlim(tint)*1e3,mvaPe?.yz.tlim(tint)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'P_e','(pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'LM','LN','MN'}',[1.02 0.9],'fontsize',fontsize);
end
if 1 % Non-gyro
  hca = irf_panel('Non-gyrotropies');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{agyro?,Dng?,Q?.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'','Non-gyrotropy'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'A_{\phi}','Dng','Q^{1/2}'}',[1.02 0.9],'fontsize',fontsize);
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.98],'color',[0 0 0],'fontsize',fontsize)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint_edr)
%irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:numel(h))
h(end).XTickLabelRotation = 0;