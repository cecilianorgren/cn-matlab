%% Plasma_density_calibration
% link to this file open('~/.d/Plasma_density_calibration.wiki')
%cd('~/data/JetBraking/20060927_1700-1800');
%tint=[toepoch([2006 9 27 17 10 0]) toepoch([2006 9 27 17 40 0])];
%tint=[toepoch([2007 9 2 14 32 0]) toepoch([2007 9 2 14 33 0])];

if 1, % read in all data
  sc_list=[3];
  for ic=sc_list,% read FGM B data and get GSM varunits
    %dobjname=irf_ssub('C?_CP_FGM_FULL',ic);eval(['caa_load ' dobjname]);varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic);
    dobjname=irf_ssub('C?_CP_FGM_FULL',ic);
    %dobjname=irf_ssub('C?_CP_FGM_SPIN',ic);
    caa_load(dobjname);
    varname=irf_ssub('B_vec_xyz_gse__C?_CP_FGM_FULL',ic);
    c_eval(['B?=getmat(' dobjname ',''' varname ''');'],ic);
    c_eval('B?=irf_abs(B?);',ic);
    c_eval('gsmB?=irf_gse2gsm(B?);',ic);
  end
  for ic=sc_list,% PEACE calculate density [cc] from PITCH_SPIN products using s/c potential correction
    c_eval('caa_load C?_CP_PEA_PITCH_SPIN_DPFlux',ic); % to speed up later
    varname=irf_ssub('Spacecraft_potential__C?_CP_EFW_L3_P',ic);
    [~,~,scpotmat]=c_caa_var_get(varname);
    scpotmat(isnan(scpotmat(:,2)),:)=[]; % remove NaN densities
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',ic);
    [var,~,varmat]=c_caa_var_get(varname);
    scpot=irf_resamp(scpotmat,varmat); % interpolate sc potential to PEACE data points
    phi=c_caa_var_get(var.DEPEND_1);
    x=getfield(c_caa_var_get(phi.DELTA_PLUS),'data');phi_dplus=x(1,:);
    x=getfield(c_caa_var_get(phi.DELTA_MINUS),'data');phi_dminus=x(1,:);
    en=c_caa_var_get(var.DEPEND_2);
    x=getfield(c_caa_var_get(en.DELTA_PLUS),'data');en_dplus=x(1,:);
    x=getfield(c_caa_var_get(en.DELTA_MINUS),'data');en_dminus=x(1,:);
    PEACE_energy_channels=en.data(1,:)+0.5*(en_dplus-en_dminus);
    PEACE_phi_min=phi.data(1,:)-phi_dminus;
    PEACE_phi_max=phi.data(1,:)+phi_dplus;
    phi_factor=repmat((cos(PEACE_phi_min*pi/180)-cos(PEACE_phi_max*pi/180))',1,length(PEACE_energy_channels));
    
    nPEACE=[varmat.t(:) varmat.t(:)*0];
    varmat.data(isnan(varmat.data))=0;
    for jj=1:size(nPEACE,1),
      satpot=(-scpot(jj,2)+2)*1.23; % assumes that probe t spacecraft potential is ~75% of spacecraft potential
      ii_energy=find(PEACE_energy_channels>satpot); % use only these energy chanels
      [en_min,ii_energy_min]=min(PEACE_energy_channels(ii_energy));
      ii_energy(ii_energy==ii_energy_min)=[]; % remove first channel after satellite potential
      [en_min,ii_energy_min]=min(PEACE_energy_channels(ii_energy));
      ii_energy(ii_energy==ii_energy_min)=[]; % remove second channel after satellite potential
      en_factor=repmat(1e-3*(en_dplus(ii_energy)+en_dminus(ii_energy)).*sqrt(1e-3*(PEACE_energy_channels(ii_energy)-satpot))./(1e-3*PEACE_energy_channels(ii_energy)),length(phi_dplus),1);
      ncoef=ones(length(phi_dplus),length(ii_energy));
      ncoef=ncoef*0.2284e-7*sqrt(1/1836)*2*pi.*phi_factor(:,ii_energy).*en_factor;
      nPEACE(jj,2)=sum(sum(shiftdim(varmat.data(jj,:,ii_energy)).*ncoef));
    end
    c_eval('nPEACE?=nPEACE;',ic);
  end
  for ic=sc_list,% PEACE calculate pressure [nPa] from PITCH_SPIN products
    energy_threshold=60; % only energy channels above this are used (to avoid photoelectrons)
    %c_eval('caa_load(C?_CP_PEA_PITCH_SPIN_DPFlux,''tint'',tint)',ic); % to speed up later
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',ic);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname);
    phi=c_caa_var_get(var.DEPEND_1);
    x=getfield(c_caa_var_get(phi.DELTA_PLUS),'data');phi_dplus=x(1,:);
    x=getfield(c_caa_var_get(phi.DELTA_MINUS),'data');phi_dminus=x(1,:);
    en=c_caa_var_get(var.DEPEND_2);
    x=getfield(c_caa_var_get(en.DELTA_PLUS),'data');en_dplus=x(1,:);
    x=getfield(c_caa_var_get(en.DELTA_MINUS),'data');en_dminus=x(1,:);
    PEACE_energy_channels=en.data(1,:)+0.5*(en_dplus-en_dminus);
    PEACE_phi_min=phi.data(1,:)-phi_dminus;
    PEACE_phi_max=phi.data(1,:)+phi_dplus;
    ii_energy=find(PEACE_energy_channels>energy_threshold); % use only these energy chanels
    Pcoef=ones(length(phi_dplus),length(ii_energy));
    phi_factor=repmat((cos(PEACE_phi_min*pi/180)-cos(PEACE_phi_max*pi/180))',1,length(ii_energy));
    en_factor=repmat(1e-3*(en_dplus(ii_energy)+en_dminus(ii_energy)).*sqrt(1e-3*PEACE_energy_channels(ii_energy)),length(phi_dplus),1);
    Pcoef=Pcoef*2/3*0.731026e-8*sqrt(1/1836)*2*pi.*phi_factor.*en_factor;
    
    P_PEACE=[varmat.t(:) varmat.t(:)*0];
    varmat.data(isnan(varmat.data))=0;
    for jj=1:size(P_PEACE,1),
      P_PEACE(jj,2)=sum(sum(shiftdim(varmat.data(jj,:,ii_energy)).*Pcoef));
    end
    c_eval('P_PEACE?=P_PEACE;',ic);
  end
  for ic=sc_list,% PEACE calculate temperature [keV] from PITCH_SPIN products (density and pressure should be calculted before)
    c_eval('T_PEACE?=irf_multiply(4.1609,P_PEACE?,1,nPEACE?,-1);',ic);
  end
  for ic=sc_list,% PEACE calculate temperature [keV] from MOMENTS
    caa_load PEA MOMENTS
    c_eval('Teperp?=getmat(C?_CP_PEA_MOMENTS,''Data_Temperature_ComponentPerpendicularToMagField__C?_CP_PEA_'');',ic);
    c_eval('Teperp?=irf_tappl(Teperp?,''*1e6*Units.kB/(Units.eV*1e3)'');',ic); % convert temperature to keV
  end
  ncal=1; % PEACE calibrate integrated density against WHISPER by this factor
  c_eval('ncal_PEACE?=[nPEACE?(:,1) nPEACE?(:,2)*ncal];',sc_list);
  for ic=sc_list,% EFW calculate satellite potential, corrected by Cully et al., 2007
    varname=irf_ssub('Spacecraft_potential__C?_CP_EFW_L3_P',ic);
    [~,~,vps]=c_caa_var_get(varname);
    c_eval('Vsat?=[vps(:,1) (-vps(:,2)+2)*1.23];',ic);
  end
end
if 1, % initialize figure
  fn=figure(61);
  h=irf_plot(7);
end
ic=3; % for which satellite to plot
if 1, % plot the figure
  if 1,   % PANEL: C?       FGM |B|
    hca=irf_panel('FGM |B|');
    varname=irf_ssub('B_mag__C?_CP_FGM_5VPS',ic);
    irf_plot(hca,varname,'nolabels');
    ylabel(hca,'|B| [nT]');
  end
  title(hca,['C' num2str(ic)]);
  if 1,   % PANEL: C?       PEACE PITCH_SPIN_DPFlux spectrogram omni
    % ic (number of spacecraft) should be defined
    hca=irf_panel('C1 PEACE energy spectra');
    dobjname=irf_ssub('C?_CP_PEA_PITCH_SPIN_DPFlux',ic);
    caa_load(dobjname,'tint',tint);
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',ic);
    varunits=eval(['getunits(' dobjname ',''' varname ''')']);
    %varunits='log_{10} dEF\newline keV/cm^2 s sr keV';
    irf_plot(hca,varname,'sum_dim1','colorbarlabel',varunits,'fitcolorbarlabel','nolabels');
    caxis(hca,[5.8 7.6]);
    irf_colormap(hca,'default');
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E [eV]');
  end
  if 1,   % ADD:   C?       EFW satellite potential (corrected)
    % change L3 to L2 to get full resolution instead of spin
    c_eval('Vsat=Vsat?;',ic);
    hold(hca,'on');
    irf_plot(hca,Vsat);
    hold(hca,'off');
    irf_zoom(hca,'y',[30 3e4]);
  end
  if 1,   % PANEL: C1..C4   PEACE temperature [keV]
    hca=irf_panel('PEACE T integrated');
    T_PEACE=irf_ssub('T_PEACE?',ic);
    irf_plot(hca,T_PEACE,'-');
    set(hca,'yscale','log');
    irf_zoom(hca,'y',[0.09 8]);
    ylabel(hca,'T_{e} [keV]');
  end
  if 1,   % ADD:   C?       PEACE MOMENTS temperature [keV]
    %
    % hca=irf_panel('PEACE T');
    Teperp=irf_ssub('Teperp?',ic);
    hold(hca,'on');
    irf_plot(hca,Teperp,'.','markersize',12,'color','r');
    hold(hca,'off');
  end
  if 1,   % PANEL: C1..C4   PEACE density [cc]
    hca=irf_panel('PEACE N');
    ncal_PEACE=irf_ssub('ncal_PEACE?',ic);
    irf_plot(hca,ncal_PEACE);
    set(hca,'yscale','log','ytick',[1e-2 1e-1]);
    irf_zoom(hca,'y',[0.006 0.5]);
    ylabel(hca,'N_{e} [cc]');
  end
  if 1,   % ADD:   C?       PEACE MOMENTS density [cc]
    %
    % hca=irf_panel('PEACE N');
    varname=irf_ssub('Data_Density__C?_CP_PEA_MOMENTS',ic);
    [~,~,nPEACE]=c_caa_var_get(varname);
    hold(hca,'on');
    irf_plot(hca,nPEACE,'.','markersize',12,'color','r');
    hold(hca,'off');
    ylabel(hca,'N [cc]');
  end
  if 1,   % PANEL: C?       EFW satellite potential
    % change L3 to L2 to get full resolution instead of spin
    hca=irf_panel('EFW satellite potential spin');
    Vps=irf_ssub('Spacecraft_potential__C?_CP_EFW_L3_P',ic);
    irf_plot(hca,Vps,'nolabels');
    ylabel(hca,'Sat pot [V]');
  end
  if 1,   % PANEL: C?       WHISPER spectrogram
    hca=irf_panel('WHISPER spectrogram natural');
    varname=irf_ssub('Electric_Spectral_Power_Density__C?_CP_WHI_NATURAL',ic);
    %[~,~,~,varunits]=c_caa_var_get(varname);
    varunits='(V/m)^2/Hz';
    % REMOVE 'fillspectgrogramgaps' flag in the next line if precise intervals of
    % WHISPER measurements are needed !!!!
    irf_plot(hca,varname,'tint',tint,'colorbarlabel',varunits,'fitcolorbarlabel','fillspectrogramgaps','nolabels');
    hold(hca,'on');
    c_eval('fpe=irf_plasma_calc(irf_resamp(B?,ncal_PEACE?),ncal_PEACE?,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
    irf_plot(hca,irf_tappl(fpe,'/1e3'),'-','linewidth',0.2,'color','k');
    c_eval('fpemom=irf_plasma_calc(irf_resamp(B?,nPEACE),nPEACE,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
    irf_plot(hca,irf_tappl(fpemom,'/1e3'),'.','linewidth',0.2,'color','r');
    hold(hca,'off');
    caxis(hca,[-16 -11]);
    set(hca,'yscale','log','ytick',[3 4 5 1e1 20 30 50 ]);
    irf_zoom(hca,'y',[2 12]);
  end
  if 1,   % PANEL: C?       STAFF spectrogram Ex and fce/flh lines
    hca=irf_panel('STAFF spectrogram Ex and fce/flh lines');
    varname=irf_ssub('EE_xxyy_isr2__C?_CP_STA_PSD',ic);
    %[~,~,~,varunits]=c_caa_var_get(varname);
    varunits='(mV/m)^2/Hz';
    irf_plot(hca,varname,'tint',tint,'colorbarlabel',varunits,'fitcolorbarlabel','comp',1,'nolabels');
    hold(hca,'on');
    c_eval('fce=irf_plasma_calc(B?,0,0,0,0,''Fce'');',ic); % calculate electron gyrofrequency
    irf_plot(hca,fce,'-','linewidth',0.2,'color','k');
    c_eval('fpe=irf_plasma_calc(irf_resamp(B?,ncal_PEACE?),ncal_PEACE?,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
    irf_plot(hca,fpe,'-','linewidth',0.2,'color','k');
    c_eval('fpemom=irf_plasma_calc(irf_resamp(B?,nPEACE),nPEACE,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
    irf_plot(hca,fpemom,'.','linewidth',0.2,'color','r');
    c_eval('flh=irf_plasma_calc(B?,1,0,0,0,''Flh'');',ic); % calculate lower hybrid frequency (underdense case in puter magnetosphere)
    irf_plot(hca,flh,'-','linewidth',0.2,'color','k');
    hold(hca,'off');
    caxis(hca,[-9 -1]);
    set(hca,'yscale','log','ytick',[1e1 1e2 1e3]);
    irf_zoom(hca,'y',[10 4000]);
  end
  
  irf_plot_axis_align
  irf_zoom(h,'x',tint);
  irf_zoom(h(1),'y');
  irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);
  irf_timeaxis(h);
end
