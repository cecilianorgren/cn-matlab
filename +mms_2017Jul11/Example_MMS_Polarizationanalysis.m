% Perform polarization analysis on burst mode electric and magnetic fields.
% Plots spectrograms, ellipticity, wave-normal angle, planarity, degree of
% polarization (DOP), phase speed, and normalized Poynting flux along B.
% Time selections should not be too long (less than 20 seconds), 
% otherwise the analysis will be very slow. 
% Written by D. B. Graham.

ic = 1; % Spacecraft number

Tint = irf.tint('2017-07-11T22:33:58.00Z/2017-07-11T22:34:08.00Z'); % EDR
Tint = irf.tint('2017-07-11T22:33:18.00Z/2017-07-11T22:33:32.00Z'); % LH
% 5 s takes 1 min

%% Load data
Tintl = Tint+[-100 100];

ic = 1;
c_eval('n = ne?.tlim(Tint);',ic)
c_eval('Rxyz = gseR?.tlim(Tint);',ic)
c_eval('Bxyz = gseB?.tlim(Tint);',ic)
c_eval('Exyz = gseE?.tlim(Tint);',ic)
c_eval('Bscm = gseB?scm.tlim(Tint);',ic)

%% Polarization analysis
Units=irf_units; % read in standard units
Me=Units.me;
e=Units.e;
epso=Units.eps0;
B_SI=Bxyz.abs.data*1e-9;
N_SI=n*1e6;

Wce = e*B_SI/Me;
ecfreq = Wce/2/pi;
ecfreq01 = ecfreq*0.1;
ecfreq05 = ecfreq*0.5;
ecfreq = irf.ts_scalar(Bxyz.time,ecfreq);
ecfreq01 = irf.ts_scalar(Bxyz.time,ecfreq01);
ecfreq05 = irf.ts_scalar(Bxyz.time,ecfreq05);

Wpe = sqrt(N_SI*e^2/Me/epso); % rad/s
fpe = Wpe/2/pi; % Hz

tic
polarization = irf_ebsp(Exyz,Bscm,Bxyz,Bxyz,Rxyz,[10 4000],'polarization','fac');
toc
%%
frequency = polarization.f;
time = polarization.t;
Bsum = polarization.bb_xxyyzzss(:,:,4);
Bperp = polarization.bb_xxyyzzss(:,:,1)+polarization.bb_xxyyzzss(:,:,2);
Bpar = polarization.bb_xxyyzzss(:,:,3);
Esum = polarization.ee_xxyyzzss(:,:,4);
Eperp = polarization.ee_xxyyzzss(:,:,1)+polarization.ee_xxyyzzss(:,:,2);
Epar = polarization.ee_xxyyzzss(:,:,3);
ellipticity = polarization.ellipticity;
dop = polarization.dop;
thetak = polarization.k_tp(:,:,1);
planarity = polarization.planarity;
pfluxz = polarization.pf_xyz(:,:,3)./sqrt(polarization.pf_xyz(:,:,1).^2+polarization.pf_xyz(:,:,2).^2+polarization.pf_xyz(:,:,3).^2);
% Calculate phase speed v_ph = E/B.
vph = sqrt(Esum./Bsum)*1e6;
vphperp = sqrt(Eperp./Bperp)*1e6;

%% Remove points with very low B amplitutes
Bsumthres = 1e-8;
removepts = find(Bsum < Bsumthres);
ellipticity(removepts) = NaN;
thetak(removepts) = NaN;
dop(removepts) = NaN;
planarity(removepts) = NaN;
pfluxz(removepts) = NaN;
vph(removepts) = NaN;
vphperp(removepts) = NaN;

%% Plot
ic = 1;
c_eval('flh = flh?;',ic)
% Define blue-red colormap
rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
bgrcmap = [rr' gg' bb'];

npanels = 11;
h=irf_plot(npanels); 
isub = 1;
isspec=zeros(npanels,1);

if 1 % B
  isub = isub + 1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
end
if 0 % J curl par perp
  isub = isub + 1;
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('12'))
  irf_plot(hca,{gseJcurl.dot(gseB1.norm),abs(gseJcurl-gseJcurl.dot(gseB1.norm)*gseB1.norm.resample(gseJcurl))},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'J_{||}','|J_{\perp}|'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 0 % Lp from Eperp
  isub = isub + 1;
  hca = irf_panel('Lp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Lp?*1e-3},''comp'');',ic)
  hca.YLabel.String = {'L_p','(km)'};  
  hca.YScale = 'log';
  hca.YLim = [1e2 5e3];
end
if 0 % Te par perp Ti/Tref
  isub = isub + 1;
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
if 0 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [100 30000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,16).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
   irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % J curl
  isub = isub + 1;
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 0 % ne
  isub = isub + 1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 0 % Ve  
  isub = isub + 1;
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % B scm
  isub = isub + 1;
  hca = irf_panel('B scm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?scm.x.tlim(tint),gseB?scm.y.tlim(tint),gseB?scm.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B_{SCM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % B sum spectrogram
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('Bsum');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=Bsum;
  specrec.f_label='';
  specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
  irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  irf_legend(hca,'(a)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-8 -4])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 0 % B perp
  isub = isub + 1;
  hca = irf_panel('B vec perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?scmperp.x.tlim(tint),gseB?scmperp.y.tlim(tint),gseB?scmperp.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B_{\perp}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % B per spectrogram
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('Bper');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=Bperp;
  specrec.f_label='';
  specrec.p_label={'log_{10}B_{\perp}^{2}','nT^2 Hz^{-1}'};
  irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  %irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  irf_plot(hca,flh,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-8 -4])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 0 % B par
  isub = isub + 1;
  hca = irf_panel('B vec par');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gseB?scmpar},''comp'');',ic)
  hca.YLabel.String = {'B_{||}','(mV/m)'};
  irf_zoom(hca,'y')
end
if 0 % B par spectrogram
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('Bpar');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=Bpar;
  specrec.f_label='';
  specrec.p_label={'log_{10}B_{||}^{2}','nT^2 Hz^{-1}'};
  irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  %irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  irf_plot(hca,flh,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-8 -4])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 1 % E tot
  isub = isub + 1;
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E sum spectrogram
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('Esum');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=Esum;
  specrec.f_label='';
  specrec.p_label={'log_{10}E^{2}','mV^2 m^{-2} Hz^{-1}'};
  irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  irf_legend(hca,'(b)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  %irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  irf_plot(hca,flh,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-6 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 0 % E perp
  isub = isub + 1;
  hca = irf_panel('E vec perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x.tlim(tint),gseE?perp.y.tlim(tint),gseE?perp.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E per spectrogram
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('Eper');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=Eperp;
  specrec.f_label='';
  specrec.p_label={'log_{10}E_{\perp}^{2}','mV^2 m^{-2} Hz^{-1}'};
  irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  %irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  irf_plot(hca,flh,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-6 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 0 % E par
  isub = isub + 1;
  hca = irf_panel('E vec par');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  irf_zoom(hca,'y')
end
if 0 % E par spectrogram
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('Epar');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=Epar;
  specrec.f_label='';
  specrec.p_label={'log_{10}E_{||}^{2}','mV^2 m^{-2} Hz^{-1}'};
  irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  %irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  irf_plot(hca,flh,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-6 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 1 % Ellipticity
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('ellipt');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=ellipticity;
  specrec.f_label='';
  specrec.p_label={'Ellipticity'};
  irf_spectrogram(hca,specrec,'lin','donotfitcolorbarlabel');
  irf_legend(hca,'(c)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-1, 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,bgrcmap);
end
if 1 % thetak, propagation angle
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('thetak');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=thetak;
  specrec.f_label='';
  specrec.p_label={'\theta_{k}'};
  irf_spectrogram(hca,specrec,'lin','donotfitcolorbarlabel');
  irf_legend(hca,'(d)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[0, 90])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 1 % degree of polarization
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('dop');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=dop;
  specrec.f_label='';
  specrec.p_label={'DOP'};
  irf_spectrogram(hca,specrec,'lin','donotfitcolorbarlabel');
  irf_legend(hca,'(e)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[0, 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 1 % planarity
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('planarity');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=planarity;
  specrec.f_label='';
  specrec.p_label={'planarity'};
  irf_spectrogram(hca,specrec,'lin','donotfitcolorbarlabel');
  irf_legend(hca,'(f)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[0, 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 1 % vph, E/B
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('vph');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=vph;
  specrec.f_label='';
  specrec.p_label={'log_{10}E/B','m s^{-1}'};
  irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  irf_legend(hca,'(g)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[5 8])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 1% poynting flux, ExB
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('poynting');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=pfluxz;
  specrec.f_label='';
  specrec.p_label={'S_{||}/|S|'};
  irf_spectrogram(hca,specrec,'lin','donotfitcolorbarlabel');
  irf_legend(hca,'(h)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[1 4])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,bgrcmap);
end

% Remove grid and set background to grey
set(h(find(isspec==1)),'xgrid','off','ygrid','off')
%set(h(find(isspec==1)),'Color',0.7*[1 1 1]);


c_eval('title(h(1),''MMS? Polarization Analysis'');',ic);
set(h,'fontsize',12);
irf_plot_axis_align(h);
irf_zoom(h,'x',Tint);

irf_zoom(h(find(isspec==0)),'y');

%% Find center of where the power peaks
Tint = irf.tint('2015-10-16T10:33:45.00Z/2015-10-16T10:33:45.70Z');

specrec=struct('t',time);
specrec.f=frequency;
specrec.p=Bsum;
specrec.f_label='';
specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};

ALL_PEAKS = nan(size(specrec.p));
for it = 1:numel(specrec.t)
  [PKS,LOCS] = findpeaks(specrec.p(it,:));
  ALL_PEAKS(it,LOCS) = PKS;
end

%% Plot

h=irf_plot(2,'newfigure'); 

hca=irf_panel('Bsum');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=Bsum;
specrec.f_label='';
specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
irf_legend(hca,'(a)',[0.99 0.98],'color','w','fontsize',12)
hold(hca,'on');
irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
hold(hca,'off');
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3]);
caxis(hca,[-8 -4])
ylabel(hca,'f (Hz)','fontsize',12);

hca=irf_panel('peaks');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=ALL_PEAKS;
specrec.f_label='';
specrec.p_label={'peaks'};
irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
irf_legend(hca,'(b)',[0.99 0.98],'color','w','fontsize',12)
hold(hca,'on');
irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
hold(hca,'off');
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3]);
caxis(hca,[-8 -4])
ylabel(hca,'f (Hz)','fontsize',12);