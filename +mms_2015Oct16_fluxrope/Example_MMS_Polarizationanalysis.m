% Perform polarization analysis on burst mode electric and magnetic fields.
% Plots spectrograms, ellipticity, wave-normal angle, planarity, degree of
% polarization (DOP), phase speed, and normalized Poynting flux along B.
% Time selections should not be too long (less than 20 seconds), 
% otherwise the analysis will be very slow. 
% Written by D. B. Graham.

ic = 1; % Spacecraft number


Tint = irf.tint('2015-10-16T10:33:43.00Z/2015-10-16T10:33:52.00Z');
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
units = irf_units; % read in standard units
Me = units.me;
e = units.e;
B_SI = Bxyz.abs.data*1e-9;
N_SI = n*1e6;

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
% Esum, Bsum is in units of E^2, B^2
vph = sqrt(Esum./Bsum)*1e6; % m/s
vphperp = sqrt(Eperp./Bperp)*1e6;

fce_mat = repmat(ecfreq.resample(irf_time(time,'epoch>epochtt')).data,1,32);
freq_mat = repmat(frequency,numel(time),1);
resonant_speed = vphperp.*(freq_mat-fce_mat)./freq_mat;
resonant_energy = resonant_speed.^2*units.me/2/units.eV;

% Remove points with very low B amplitutes
Bsumthres = 1e-7;
removepts = find(Bsum < Bsumthres);
ellipticity(removepts) = NaN;
thetak(removepts) = NaN;
dop(removepts) = NaN;
planarity(removepts) = NaN;
pfluxz(removepts) = NaN;
vph(removepts) = NaN;
vphperp(removepts) = NaN;
resonant_energy(removepts) = NaN;
% Remove points outside given vph interval (to isolate waves)
vphint = [1e6 1e7];
resonant_energy(find(vph > vphint(2))) = NaN;
resonant_energy(find(vph < vphint(1))) = NaN;
resonant_energy(find(dop < 0.8)) = NaN;
resonant_energy(find(resonant_energy > 5e3)) = NaN;

%% Find resonant energy using the spectra
resonant_energy_limits = nan(numel(time),2);
resonant_energy_limits_freq = nan(numel(time),2);
resonant_energy_mid = nan(numel(time),1);
resonant_energy_mid_freq = nan(numel(time),1);
for it = 1:numel(time)
  if any(not(isnan(resonant_energy(it,:))))
    [resonant_energy_min_tmp,ind_min] = min(resonant_energy(it,:));
    [resonant_energy_max_tmp,ind_max] = max(resonant_energy(it,:));
    resonant_energy_limits(it,:) = [resonant_energy_min_tmp resonant_energy_max_tmp];    
    resonant_energy_limits_freq(it,:) = frequency([ind_min,ind_max]);
    
    if sum(not(isnan(resonant_energy(it,:))))>2
      %[resonant_energy_mid_tmp,ind_mid] = find(resonant_energy(it,:)==median(resonant_energy(it,not(isnan(resonant_energy(it,:)))),'omitnan'));
      mean_tmp = mean(resonant_energy(it,not(isnan(resonant_energy(it,:)))),'omitnan');
      median_tmp = median(resonant_energy(it,not(isnan(resonant_energy(it,:)))),'omitnan');
      
      ind_mid = find(abs(resonant_energy(it,:)-mean_tmp) == min(abs(resonant_energy(it,:)-mean_tmp)));
      %[resonant_energy_mid_tmp,ind_mid] = find(resonant_energy(it,:)==median(resonant_energy(it,not(isnan(resonant_energy(it,:)))),'omitnan'));
      resonant_energy_mid(it,:) = resonant_energy(it,ind_mid);
      resonant_energy_mid_freq(it,:) = frequency(ind_mid);    
      
      if ind_min>ind_max
       % break
      end
    end    
  end
end
ts_res_energy_mid = irf.ts_scalar(irf_time(time,'epoch>epochtt'),resonant_energy_mid);
ts_res_energy_mid_freq = irf.ts_scalar(irf_time(time,'epoch>epochtt'),resonant_energy_mid_freq);
ts_res_energy_lim = irf.ts_scalar(irf_time(time,'epoch>epochtt'),resonant_energy_limits);
ts_res_energy_lim_freq = irf.ts_scalar(irf_time(time,'epoch>epochtt'),resonant_energy_limits_freq);

% Find resonant energy using the a given average vph, it seems quite constant
% vph_choice = 10^6.5; % m/s
% frequency-ecfreq.resample(irf_time(time,'epoch>epochtt'))
% resonant_speed = vphperp.*(freq_mat-fce_mat)./freq_mat;
% resonant_energy = resonant_speed.^2*units.me/2/units.eV;

%% Plot, more adaptive
ic = 1;
figure(30)
Tint = irf.tint('2015-10-16T10:33:43.00Z/2015-10-16T10:33:52.00Z');
tint = Tint;
tintZoom = tint;

colors = mms_colors('matlab');

% Define blue-red colormap
rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
bgrcmap = [rr' gg' bb'];

npanels = 6;
h=irf_plot(npanels); 
isub = 1;
isspec=zeros(npanels,1);
zoomE = [];

if 1 % B
  isub = isub + 1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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
if 1 % eDist omni 64 apar
  hca = irf_panel('e DEF 3'); zoomE = [zoomE isub]; isub = isub + 1; 
  pas = [150 180];
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('ePitch?lim = ePDist?.tlim(tintZoom+[-2 2]).pitchangles(dmpaB?,pas);',ic)      
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % add resonant energy
    hold(hca,'on')
    irf_plot(hca,{ts_res_energy_mid},'color',colors(3,:))
    hold(hca,'off')
  end
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  if 0 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')        
    set(hca,'ColorOrder',mms_colors('1'))
    c_eval('lineWperp = irf_plot(hca,{Wpar*Efactorpar},''comp'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    %irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactorpar)},[0.77 0.2])
    hold(hca,'off')
    %irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 perp
  hca = irf_panel('e DEF 2'); zoomE = [zoomE isub]; isub = isub + 1;
  pas = [75 105];
  c_eval('ePitch?lim = ePDist?.tlim(tintZoom+[-2 2]).pitchangles(dmpaB?,pas);',ic)    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  if 0 % plot perpendicular energy baed on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('1'))
    c_eval('lineWperp = irf_plot(hca,{Wperp},''comp'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 par
  hca = irf_panel('e DEF 1'); zoomE = [zoomE isub]; isub = isub + 1;
  pas = [0 30]; 
  c_eval('ePitch?lim = ePDist?.tlim(tintZoom+[-2 2]).pitchangles(dmpaB?,pas);',ic)     
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  if 0 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')        
    set(hca,'ColorOrder',mms_colors('1'))
    c_eval('lineWperp = irf_plot(hca,{Wpar*Efactorpar},''comp'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    %irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactorpar)},[0.77 0.2])
    hold(hca,'off')
    %irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % B scm
  isub = isub + 1;
  hca = irf_panel('B scm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?scm.x.tlim(tint),gseB?scm.y.tlim(tint),gseB?scm.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B_{SCM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % B sum spectrogram
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
  irf_legend(hca,'(b)',[0.99 0.98],'color','w','fontsize',12)
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
  %%
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('Bper');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=Bperp;
  specrec.f_label='';
  specrec.p_label={'log_{10}B_{\perp}^{2}','nT^2 Hz^{-1}'};
  irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  irf_legend(hca,'(b)',[0.99 0.98],'color','w','fontsize',12)
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
if 0 % E tot
  isub = isub + 1;
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E sum spectrogram
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
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
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
  irf_legend(hca,'(b)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
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
  irf_legend(hca,'(b)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.5,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.5,'color','w')
  irf_plot(hca,fpe,'linewidth',1.5,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-6 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 0 % Ellipticity
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
if 0 % thetak, propagation angle
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
if 0 % degree of polarization
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
if 0 % planarity
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
if 1 % resonant energy
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('resonant energy');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=resonant_energy;
  specrec.f_label='';
  specrec.p_label={'log_{10}E_{res}',''};
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
  caxis(hca,[1 4])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end
if 0 % poynting flux, ExB
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
  caxis(hca,[-1 1])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,bgrcmap);
end

% Remove grid and set background to grey
set(h(find(isspec==1)),'xgrid','off','ygrid','off')
set(h(find(isspec==1)),'Color',0.7*[1 1 1]);


c_eval('title(h(1),''MMS? Polarization Analysis'');',ic);
irf_plot_axis_align(h);
irf_zoom(h,'x',Tint);
set(h,'fontsize',12);

irf_zoom(h(find(isspec==0)),'y');
irf_zoom(h(zoomE),'y',[10 2000]);

%% Plot, old

h=irf_plot(8,'newfigure'); 
%h=irf_figure(540+ic,8);
xSize=750; ySize=800;
set(gcf,'Position',[10 10 xSize ySize]);

xwidth = 0.88;
ywidth = 0.115;
set(h(1),'position',[0.08 0.975-ywidth xwidth ywidth]);
set(h(2),'position',[0.08 0.975-2*ywidth xwidth ywidth]);
set(h(3),'position',[0.08 0.975-3*ywidth xwidth ywidth]);
set(h(4),'position',[0.08 0.975-4*ywidth xwidth ywidth]);
set(h(5),'position',[0.08 0.975-5*ywidth xwidth ywidth]);
set(h(6),'position',[0.08 0.975-6*ywidth xwidth ywidth]);
set(h(7),'position',[0.08 0.975-7*ywidth xwidth ywidth]);
set(h(8),'position',[0.08 0.975-8*ywidth xwidth ywidth]);

h(1)=irf_panel('Bsum');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=Bsum;
specrec.f_label='';
specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
irf_spectrogram(h(1),specrec,'log','donotfitcolorbarlabel');
irf_legend(h(1),'(a)',[0.99 0.98],'color','w','fontsize',12)
hold(h(1),'on');
irf_plot(h(1),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(1),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(1),ecfreq01,'linewidth',1.5,'color','w')
hold(h(1),'off');
set(h(1),'yscale','log');
set(h(1),'ytick',[1e1 1e2 1e3]);
caxis(h(1),[-8 -4])
ylabel(h(1),'f (Hz)','fontsize',12);

h(2)=irf_panel('Esum');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=Esum;
specrec.f_label='';
specrec.p_label={'log_{10}E^{2}','mV^2 m^{-2} Hz^{-1}'};
irf_spectrogram(h(2),specrec,'log','donotfitcolorbarlabel');
irf_legend(h(2),'(b)',[0.99 0.98],'color','w','fontsize',12)
hold(h(2),'on');
irf_plot(h(2),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(2),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(2),ecfreq01,'linewidth',1.5,'color','w')
hold(h(2),'off');
set(h(2),'yscale','log');
set(h(2),'ytick',[1e1 1e2 1e3]);
caxis(h(2),[-6 1])
ylabel(h(2),'f (Hz)','fontsize',12);

h(3)=irf_panel('ellipt');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=ellipticity;
specrec.f_label='';
specrec.p_label={'Ellipticity'};
irf_spectrogram(h(3),specrec,'lin','donotfitcolorbarlabel');
irf_legend(h(3),'(c)',[0.99 0.98],'color','w','fontsize',12)
hold(h(3),'on');
irf_plot(h(3),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(3),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(3),ecfreq01,'linewidth',1.5,'color','w')
hold(h(3),'off');
set(h(3),'yscale','log');
set(h(3),'ytick',[1e1 1e2 1e3]);
caxis(h(3),[-1, 1])
ylabel(h(3),'f (Hz)','fontsize',12);

h(4)=irf_panel('thetak');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=thetak;
specrec.f_label='';
specrec.p_label={'\theta_{k}'};
irf_spectrogram(h(4),specrec,'lin','donotfitcolorbarlabel');
irf_legend(h(4),'(d)',[0.99 0.98],'color','w','fontsize',12)
hold(h(4),'on');
irf_plot(h(4),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(4),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(4),ecfreq01,'linewidth',1.5,'color','w')
hold(h(4),'off');
set(h(4),'yscale','log');
set(h(4),'ytick',[1e1 1e2 1e3]);
caxis(h(4),[0, 90])
ylabel(h(4),'f (Hz)','fontsize',12);

h(5)=irf_panel('dop');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=dop;
specrec.f_label='';
specrec.p_label={'DOP'};
irf_spectrogram(h(5),specrec,'lin','donotfitcolorbarlabel');
irf_legend(h(5),'(e)',[0.99 0.98],'color','w','fontsize',12)
hold(h(5),'on');
irf_plot(h(5),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(5),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(5),ecfreq01,'linewidth',1.5,'color','w')
hold(h(5),'off');
set(h(5),'yscale','log');
set(h(5),'ytick',[1e1 1e2 1e3]);
caxis(h(5),[0, 1])
ylabel(h(5),'f (Hz)','fontsize',12);

h(6)=irf_panel('planarity');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=planarity;
specrec.f_label='';
specrec.p_label={'planarity'};
irf_spectrogram(h(6),specrec,'lin','donotfitcolorbarlabel');
irf_legend(h(6),'(f)',[0.99 0.98],'color','w','fontsize',12)
hold(h(6),'on');
irf_plot(h(6),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(6),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(6),ecfreq01,'linewidth',1.5,'color','w')
hold(h(6),'off');
set(h(6),'yscale','log');
set(h(6),'ytick',[1e1 1e2 1e3]);
caxis(h(6),[0, 1])
ylabel(h(6),'f (Hz)','fontsize',12);
  
h(7)=irf_panel('vph');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=vph;
specrec.f_label='';
specrec.p_label={'log_{10}E/B','m s^{-1}'};
irf_spectrogram(h(7),specrec,'log','donotfitcolorbarlabel');
irf_legend(h(7),'(g)',[0.99 0.98],'color','w','fontsize',12)
hold(h(7),'on');
irf_plot(h(7),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(7),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(7),ecfreq01,'linewidth',1.5,'color','w')
hold(h(7),'off');
set(h(7),'yscale','log');
set(h(7),'ytick',[1e1 1e2 1e3]);
caxis(h(7),[5 8])
ylabel(h(7),'f (Hz)','fontsize',12);

h(8)=irf_panel('poynting');
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=pfluxz;
specrec.f_label='';
specrec.p_label={'S_{||}/|S|'};
irf_spectrogram(h(8),specrec,'lin','donotfitcolorbarlabel');
irf_legend(h(8),'(h)',[0.99 0.98],'color','w','fontsize',12)
hold(h(8),'on');
irf_plot(h(8),ecfreq,'linewidth',1.5,'color','w')
irf_plot(h(8),ecfreq05,'linewidth',1.5,'color','w')
irf_plot(h(8),ecfreq01,'linewidth',1.5,'color','w')
hold(h(8),'off');
set(h(8),'yscale','log');
set(h(8),'ytick',[1e1 1e2 1e3]);
caxis(h(8),[-1 1])
ylabel(h(8),'f (Hz)','fontsize',12);

% Remove grid and set background to grey
set(h(1:8),'xgrid','off','ygrid','off')
set(h(1:8),'Color',0.7*[1 1 1]);

% Define blue-red colormap
rr = interp1([1 64 128 192 256],[0.0  0.5 0.75 1.0 0.75],1:256);
gg = interp1([1 64 128 192 256],[0.0  0.5 0.75 0.5 0.00],1:256);
bb = interp1([1 64 128 192 256],[0.75 1.0 0.75 0.5 0.00],1:256);
bgrcmap = [rr' gg' bb'];

colormap(h(1),'jet');
colormap(h(2),'jet');
colormap(h(3),bgrcmap);
colormap(h(4),'jet');
colormap(h(5),'jet');
colormap(h(6),'jet');
colormap(h(7),'jet');
colormap(h(8),bgrcmap);

c_eval('title(h(1),''MMS? Polarization Analysis'');',ic);
irf_plot_axis_align(h(1:8));
irf_zoom(h(1:8),'x',Tint);
set(h(1:8),'fontsize',12);

%% Find center of where the power peaks
Tint = irf.tint('2015-10-16T10:33:45.00Z/2015-10-16T10:33:45.70Z');


specrec=struct('t',time);
specrec.f=frequency;
specrec.p=Bsum;
specrec.f_label='';
specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};

% remove points to isolate whistler waves
Bsumthres = 1e-7;
removepts = find(Bsum < Bsumthres);
specrec.p(removepts) = NaN;
ellipticitythres = 0.8;
removepts = find(ellipticity < ellipticitythres);
specrec.p(removepts) = NaN;

ALL_PEAKS = zeros(size(specrec.p));
resonant_energy_peaks = nan(numel(specrec.t),1);
for it = 1:numel(specrec.t)
  [PKS,LOCS] = findpeaks(specrec.p(it,:));    
  ALL_PEAKS(it,LOCS) = PKS;
  [val,loc] = find(PKS == max(PKS));
  if not(isempty(val))
    resonant_energy_peaks(it,1) = val;
  end
end


%% Find ellipticity above certain value
specrec=struct('t',time);
specrec.f=frequency;
specrec.p=ellipticity;
specrec.f_label='';
specrec.p_label={'Ellipticity'};

ALL_PEAKS = nan(size(specrec.p));
for it = 1:numel(specrec.t)
  [PKS,LOCS] = find(specrec.p(it,:)>0.8);
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