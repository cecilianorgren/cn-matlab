% Perform polarization analysis on burst mode electric and magnetic fields.
% Plots spectrograms, ellipticity, wave-normal angle, planarity, degree of
% polarization (DOP), phase speed, and normalized Poynting flux along B.
% Time selections should not be too long (less than 20 seconds), 
% otherwise the analysis will be very slow. 
% Written by D. B. Graham.

ic = 1; % Spacecraft number

Tint = irf.tint('2017-07-11T22:33:58.00Z/2017-07-11T22:34:08.00Z'); %20151112071854
% 5 s takes 1 min

%% Load data
Tintl = Tint+[-100 100];

Rxyz = gseR1.tlim(Tint);
Bxyz = gseB1.tlim(Tint);
Exyz = gseE1.tlim(Tint);
Bscm = gseB1scm.tlim(Tint);

%% Polarization analysis
Units=irf_units; % read in standard units
Me=Units.me;
e=Units.e;
B_SI=Bxyz.abs.data*1e-9;
Wce = e*B_SI/Me;
ecfreq = Wce/2/pi;
ecfreq01 = ecfreq*0.1;
ecfreq05 = ecfreq*0.5;
ecfreq = irf.ts_scalar(Bxyz.time,ecfreq);
ecfreq01 = irf.ts_scalar(Bxyz.time,ecfreq01);
ecfreq05 = irf.ts_scalar(Bxyz.time,ecfreq05);

tic
polarization = irf_ebsp(Exyz,Bscm,Bxyz,Bxyz,Rxyz,[10 4000],'polarization','fac');
toc
%%
frequency = polarization.f;
time = polarization.t;
Bsum = polarization.bb_xxyyzzss(:,:,4);
Bperp = polarization.bb_xxyyzzss(:,:,1)+polarization.bb_xxyyzzss(:,:,2);
Esum = polarization.ee_xxyyzzss(:,:,4);
Eperp = polarization.ee_xxyyzzss(:,:,1)+polarization.ee_xxyyzzss(:,:,2);
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
% Define blue-red colormap
rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
bgrcmap = [rr' gg' bb'];

npanels = 12;
h=irf_plot(npanels,'newfigure'); 


if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 1 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % B spectrogram
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
  colormap(hca,'jet');
end
if 1 % E spectrogram
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
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-6 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 1 % Ellipticity
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
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-1, 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,bgrcmap);
end
if 1 % thetak, propagation angle
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
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[0, 90])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 1 % degree of polarization
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
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[0, 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 1 % planarity
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
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[0, 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 1 % vph, E/B
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
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[6 9])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 1 % poynting flux, ExB
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
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-1 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,bgrcmap);
end

% Remove grid and set background to grey
set(h(5:npanels),'xgrid','off','ygrid','off')
%set(h(1:8),'Color',0.7*[1 1 1]);




colormap(hca,'jet');
colormap(hca,bgrcmap);
colormap(hca,'jet');
colormap(hca,'jet');
colormap(hca,'jet');
colormap(hca,'jet');
colormap(hca,bgrcmap);

c_eval('title(h(1),''MMS? Polarization Analysis'');',ic);
irf_plot_axis_align(h);
irf_zoom(h,'x',Tint);
set(h,'fontsize',12);

irf_zoom(h(1:4),'y');

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