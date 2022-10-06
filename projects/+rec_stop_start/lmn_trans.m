% Set up new coordinate system
tint_mva = irf.tint('2017-07-25T20:14:08.398745849Z/2017-07-25T21:57:25.093696533Z'); % early
%tint_mva = irf.tint('2017-07-25T21:36:16.246635253Z/2017-07-25T23:52:59.490427001Z'); % later

[out,l,v]=irf_minvar(gseB1_srvy.tlim(tint_mva));

Radjust = [-1 0 0; 0 -1 0;0 0 1]; % early
%Radjust = [-1 0 0; 0 1 0;0 0 -1]; % later
%Radjust = [1 0 0; 0 1 0; 0 0 1];
R = v*Radjust';

%% USE THIS ONE FOR LARGE SCALE
r1 = [1 -0.3 0]; r1 = r1/norm(r1);
r3 = cross(cross(r1,[0 0 1]),r1);
%r3 = [0 0 1];
r2 = cross(r3,r1);
R = [r1;r2;r3];

%% Density cavity coordinates
r1 = -[-0.85 -0.52 -0.10]; r1 = r1/norm(r1);
r3 = cross(cross(r1,-[-0.21 -0.12 -0.97]),r1);
r2 = cross(r3,r1);
R = [r1;r2;r3];

%% different event
r1 = [1 0 0.5]; r1 = r1/norm(r1);
% r2 = cross(cross(r1,[0 1 0]),r1);
% r3 = cross(r1,r2);
r3 = cross(cross(r1,[-0.56 -0.46 0.68]),r1);
r2 = cross(r3,r1);

R = [r1;r2;r3];

%%
L = R(1,:);
M = R(2,:);
N = R(3,:);

% Brst
c_eval('lmnB? = gseB?*R'';',ic)
c_eval('lmnE? = gseE?*R'';',ic)
c_eval('lmnVi? = gseVi?*R'';',ic)
c_eval('lmnVe? = gseVe?*R'';',ic)
c_eval('lmnVExB? = gseVExB?*R'';',ic)
c_eval('lmnJxB? = gseJxB?*R'';',ic)

c_eval('lmnVi? = gseVi?*R'';',ic)
c_eval('lmnJ? = gseJ?*R'';',ic)
c_eval('lmnJe? = gseJe?*R'';',ic)
c_eval('lmnJi? = gseJi?*R'';',ic)


c_eval('lmnVExB? = gseVExB?*R'';',ic)

c_eval('lmnE?perp = gseE?perp*R'';',ic)
c_eval('lmnVi?perp = gseVi?perp*R'';',ic)
c_eval('lmnVi?par = gseVi?par'';',ic)

c_eval('lmnVe?perp = gseVe?perp*R'';',ic)
c_eval('lmnVe?par = gseVe?par'';',ic)

c_eval('lmnVe?perp = gseVe?perp*R'';',ic)
c_eval('lmnVi?perp = gseVi?perp*R'';',ic)


c_eval('lmnVi? = gseVi?*R'';',ic)


c_eval('lmnJxBne?_mVm = gseJxBne?_mVm*R'';',ic)
lmnJcurl = gseJcurl*R';
lmnJxB = JxB*R';
lmnJav = gseJav*R';
lmnJiav = gseJiav*R';
lmnJeav = gseJeav*R';
lmnBav = gseBav*R';

% Srvy/fast
c_eval('lmnB?_srvy = gseB?_srvy*R'';',ic)
c_eval('lmnE?_fast = gseE?_fast*R'';',ic)
c_eval('lmnVi?_fast = gseVi?_fast*R'';',ic)
c_eval('lmnVe?_fast = gseVe?_fast*R'';',ic)
c_eval('lmnVExB?_srvy = gseVExB?_srvy*R'';',ic)

%% Plot
npanels = 7;
h = irf_plot(npanels);
doFilt = 0;
ffilt = 1;

if 1 % B lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B lmn');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{lmnB?.x,lmnB?.y,lmnB?.z,lmnB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B^{LMN}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'  [L; M; N] = ',...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(1,1),R(1,2),R(1,3)),...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(2,1),R(2,2),R(2,3)),...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(3,1),R(3,2),R(3,3))}',[1.002 0.9],'fontsize',12,'color','k');  
end
if 1 % E lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E lmn');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  if doFilt
    c_eval('irf_plot(hca,{lmnE?.filt(0,ffilt,[],5).x,lmnE?.filt(0,ffilt,[],5).y,lmnE?.filt(0,ffilt,[],5).z},''comp'');',ic)    
    irf_legend(hca,sprintf('f < %g Hz',ffilt),[0.02 0.07])
  else
    c_eval('irf_plot(hca,{lmnE?.x,lmnE?.y,lmnE?.z},''comp'');',ic)    
  end
  hca.YLabel.String = {'E^{LMN}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end 
if 0 % Vi gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?_fast.x,gseVi?_fast.y,gseVi?_fast.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i^{GSE}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{lmnVi?.x,lmnVi?.y,lmnVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i^{LMN}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{lmnVi?perp.x,lmnVi?perp.y,lmnVi?perp.z,lmnVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,L}','v_{\perp,M}','v_{\perp,N}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?_fast.x,gseVe?_fast.y,gseVe?_fast.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e^{GSE}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('vplot = lmnVe?_fast;',ic)
  vplot.data(vplot.abs.data>1e4,:) = NaN;
  c_eval('irf_plot(hca,{vplot.x,vplot.y,vplot.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e^{LMN}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % VExB gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?_srvy.x,gseVExB?_srvy.y,gseVExB?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}^{GSE}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % VExB lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{lmnVExB?.x,lmnVExB?.y,lmnVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}^{LMN}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % VExB gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?_srvy.x,gseVExB?_srvy.y,gseVExB?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % n log scale
  hca = irf_panel('n log scale');
  set(hca,'ColorOrder',mms_colors('1234'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.8])
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YScale = 'log';
  hca.YLim(1) = 1e-5;
end
if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_fast.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 1 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');
    hexb_o = irf_plot(hca,Eexb*16,'color',0.5*[1 1 1],'linestyle','-');
    hold(hca,'off')
    irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?_fast_par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end

if 1 % Lp
  hca = irf_panel('Lp');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{Lp},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  %irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'L_p','(...)'};
  hca.YLim = [0 10000];
end

irf_zoom(h,'x',tint)
irf_plot_axis_align