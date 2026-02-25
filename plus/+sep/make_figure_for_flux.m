%% Load datastore
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');

events = 1:20;
nevents = numel(events);
saveAccPotPath = ['/Users/' localuser '/MATLAB/cn-matlab/+sep/acc_potential/'];
printAccPotPath = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/acceleration_potential/'];

for ievent = 1:nevents
  event = events(ievent); % set event id
  sep.get_tints; % get tints for that event id
 
  %% Load data  
  ic = 1:4;
  units = irf_units;
  % Core data
  c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
  if isempty(dmpaB4); ic = 1:3; end
  c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
  c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
  c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)  
    
  % Auxillary data, for context
  c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
  c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
  c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
  c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
  c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);

  % Auxillary quantities, for context
  c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
  c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)  
  c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
  c_eval('beta?e = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
  c_eval('beta?i = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
  c_eval('beta? = beta?i + beta?e.resample(beta?i);',ic)
  c_eval('vte?par = (2*units.eV*facTe?.xx/units.me).^.5*1e-3; vte?par.name = ''vtepar''; vte?par.units = ''km/s'';',ic)
  c_eval('vte?perp = (2*units.eV*0.5*(facTe?.yy+facTe?.zz)/units.me).^.5*1e-3; vte?perp.name = ''vteperp''; vte?perp.units = ''km/s'';',ic)
  


  %% Pick out reconnection quadrant: tail, Earth to the left
  quadrant = ones(2,2);  
  c_eval('vix? = mean(gseVi?.x.data,1);',ic)
  if ic(end) == 3
    vix = mean([vix1, vix2, vix3]);
  else
    vix = mean([vix1, vix2, vix3, vix4]);
  end
  if vix > 0; quadrant(:,2) = 0; else quadrant(:,1) = 0; end
  
  c_eval('bix? = mean(dmpaB?.tlim(tint_sep).x.data,1);',ic)
  if ic(end) == 3
    bix = mean([bix1, bix2, bix3]);
  else
    bix = mean([bix1, bix2, bix3, bix4]);
  end
  if bix > 0; quadrant(2,:) = 0; else  quadrant(1,:) = 0; end
  
  %% Pick out parameters in the lobe, sheet, and separatrix intervals
  ic_orig = ic;
  ic = 1;
  c_eval('B_lobe = mean(gseB?.abs.tlim(tint_lobe).data,1);',ic)

  c_eval('n_lobe = mean(ne?.tlim(tint_lobe).data,1);',ic)
  c_eval('n_sheet = mean(ne?.tlim(tint_sheet).data,1);',ic)
  c_eval('n_sep = mean(ne?.tlim(tint_sep).data,1);',ic)
  c_eval('n_sep_min = min(ne?.tlim(tint_sep).data,1);',ic)
  c_eval('n_sep_data = ne?.tlim(tint_sep);',ic)

  c_eval('be_lobe = mean(beta?e.tlim(tint_lobe).data,1);',ic)
  c_eval('be_sheet = mean(beta?e.tlim(tint_sheet).data,1);',ic)
  c_eval('be_sep = mean(beta?e.tlim(tint_sep).data,1);',ic)

  c_eval('Tepar_lobe = mean(facTe?.tlim(tint_lobe).data(:,1,1),1);',ic)
  c_eval('Tepar_sheet = mean(facTe?.tlim(tint_sheet).data(:,1,1),1);',ic)
  c_eval('Tepar_sep = mean(facTe?.tlim(tint_sep).data(:,1,1),1);',ic)

  c_eval('Teperp_lobe  = 0.5*mean(facTe?.tlim(tint_lobe).data(:,2,2)+facTe?.tlim(tint_lobe).data(:,3,3),1);',ic)
  c_eval('Teperp_sheet = 0.5*mean(facTe?.tlim(tint_sheet).data(:,2,2)+facTe?.tlim(tint_sheet).data(:,3,3),1);',ic)
  c_eval('Teperp_sep   = 0.5*mean(facTe?.tlim(tint_sep).data(:,2,2)+facTe?.tlim(tint_sep).data(:,3,3),1);',ic)

  wce0 = units.e*B_lobe*1e-9/units.me; % fce0*2pi
  wpe0 = sqrt(n_sheet*1e6*units.e^2/units.eps0/units.me); % fpe0*2*pi

  info_str = {sprintf('B_{lobe}=%.0f nT',B_lobe);...
              sprintf('n_{lobe}=%.3f cc',n_lobe);...
              sprintf('n_{sheet}=%.3f cc',n_sheet);...
              sprintf('n_{sep}=%.3f cc',n_sep);...
              sprintf('T_{e,par,lobe}=%.0f eV',Tepar_lobe);...
              sprintf('T_{e,par,sheet}=%.0f eV',Tepar_sheet);...
              sprintf('T_{e,par,sep}=%.0f eV',Tepar_sep);...
              sprintf('T_{e,perp,lobe}=%.0f eV',Teperp_lobe);...
              sprintf('T_{e,perp,sheet}=%.0f eV',Teperp_sheet);...
              sprintf('T_{e,perp,sep}=%.0f eV',Teperp_sep);...
              sprintf('f_{ce,0}=%.0f Hz',wce0/2/pi);...
              sprintf('f_{pe,0}=%.0f Hz',wpe0/2/pi);...
              sprintf('f_{pe,0}/f_{ce,0}=%.1f',wpe0/wce0);...
              sprintf('beta_{e,lobe}=%.3f',be_lobe);...
              sprintf('beta_{e,sheet}=%.3f',be_sheet);...
              sprintf('beta_{e,sep}=%.3f',be_sep);...
              };
  colors = mms_colors('xzy');
  info_color = [0 0 0; colors; colors; colors; 0 0 0; 0 0 0; 0 0 0; colors];

  tint_day_utc = tint.utc; tint_day_utc = tint_day_utc(1,1:10);
  tint_lobe_utc = tint_lobe.utc; tint_lobe_utc = tint_lobe_utc(:,12:23);
  tint_sheet_utc = tint_sheet.utc; tint_sheet_utc = tint_sheet_utc(:,12:23);
  tint_sep_utc = tint_sep.utc; tint_sep_utc = tint_sep_utc(:,12:23);
  tint_phi_utc = tint_phi.utc; tint_phi_utc = tint_phi_utc(:,12:23);

  str_print=sprintf('sep_ov_%s_%s_%s',tint_day_utc,tint_sep_utc(1,:),tint_sep_utc(2,:));
  str_print(strfind(str_print,':'))=[];

  ic = ic_orig;
  
  %% Remove photoelectron and secondary electron background, for tint_phi
  %tint_fred = tint_phi;%tint_phi;
  tint_phi_ = tint_fred;
  c_eval('eDist?_orig = ePDist?.tlim(tint_phi_);',ic)  
  % Remove background, try first with these default values, mostly for
  % cosmetics anyways
  nSecondary = 5;
  nPhoto = 1;
  c_eval('eDist?_nobg = mms.remove_edist_background(eDist?_orig,''nSecondary'',nSecondary,''Nphotoe_art'',nPhoto);',ic)

  %% Compute reduced electron distribution
  eint = [00 40000];
  lowerelim = 000; % not used when removing background
      
  vgmax = 70000;
  vg = -vgmax:1000:vgmax;
  vg(abs(vg)>70000) = [];

  nMC = 200;
  
  c_eval('tic; ef1D?_orig = eDist?_orig.reduce(''1D'',dmpaB?.resample(eDist?_orig).norm,''scpot'',scPot?.resample(eDist?_orig),''lowerelim'',lowerelim,''nMC'',nMC,''vg'',vg); toc',ic) % reduced distribution along B
  c_eval('tic; ef1D?_nobg = eDist?_nobg.reduce(''1D'',dmpaB?.resample(eDist?_nobg).norm,''scpot'',scPot?.resample(eDist?_nobg),''lowerelim'',lowerelim,''nMC'',nMC,''vg'',vg); toc',ic) % reduced distribution along B

  %% Get acceleration potential
  c_eval('f_lobe = ef1D?_nobg.tlim(tint_lobe);',1)
  c_eval('f_sheet = ef1D?_nobg.tlim(tint_sheet);',1)
  f_lobe = max(mean(f_lobe.data,1));
  f_sheet = max(mean(f_sheet.data,1));
  f_prom = 0.05*f_lobe*n_sep/n_lobe;
  
  eint = [000 40000];
  if any(quadrant([1 4])) % top left or bottom right, inflow is antiparallel to B
    vint = [-100000 0];
  else  % top right or bottom left, inflow is parallel to B
    vint = [0 100000];
  end
  relativeminpeakprominence = 0.15;
  minpeakprominence = f_prom;
  resample_timestep = 3;
  % Need to make for-loop for printing, since figure is done within script
  for iic = ic    
    c_eval('tsAccPot?_orig_rel = find_acc_pot(ef1D?_orig,''eint'',eint,''vint'',vint,''relativeminpeakprominence'',relativeminpeakprominence,''resample'',resample_timestep);',iic)    
    hcf = gcf;
    hcf.Position = [0 100 700 800];
    h = findobj(hcf.Children,'type','axes');
    c_eval('irf_pl_mark(h(?),tint_phi,''r'');',1:numel(h))
    cn.print(sprintf('event%g_tsAccPot%g_orig_rel',event,iic),'path',printAccPotPath)
    
    c_eval('tsAccPot?_orig_abs = find_acc_pot(ef1D?_orig,''eint'',eint,''vint'',vint,''minpeakprominence'',minpeakprominence,''resample'',resample_timestep);',iic)    
    hcf = gcf;
    hcf.Position = [0 100 700 800];
    cn.print(sprintf('event%g_tsAccPot%g_orig_abs',event,iic),'path',printAccPotPath)
    
    c_eval('tsAccPot?_nobg_rel = find_acc_pot(ef1D?_nobg,''eint'',eint,''vint'',vint,''relativeminpeakprominence'',relativeminpeakprominence,''resample'',resample_timestep);',iic)    
    hcf = gcf;   
    hcf.Position = [0 100 700 800];
    cn.print(sprintf('event%g_tsAccPot%g_nobg_rel',event,iic),'path',printAccPotPath)
    
    
    c_eval('tsAccPot?_nobg_abs = find_acc_pot(ef1D?_nobg,''eint'',eint,''vint'',vint,''minpeakprominence'',minpeakprominence,''resample'',resample_timestep);',iic)    
    hcf = gcf;    
    hcf.Position = [0 100 700 800];
    cn.print(sprintf('event%g_tsAccPot%g_nobg_abs',event,iic),'path',printAccPotPath)
    
  end
  clear tsAccPot_orig_rel tsAccPot_orig_abs tsAccPot_nobg_rel tsAccPot_nobg_abs
  c_eval('tsAccPot_orig_rel{?} = tsAccPot?_orig_rel;',ic)
  c_eval('tsAccPot_orig_abs{?} = tsAccPot?_orig_abs;',ic)
  c_eval('tsAccPot_nobg_rel{?} = tsAccPot?_nobg_rel;',ic)
  c_eval('tsAccPot_nobg_abs{?} = tsAccPot?_nobg_abs;',ic)
          
  % From timeseries, pick out maximal value
  % Can do this when reloading data
  c_eval('acc_pot?_orig_abs = max(tsAccPot?_orig_abs.data);',ic)
  c_eval('acc_pot?_nobg_abs = max(tsAccPot?_nobg_abs.data);',ic)
  c_eval('acc_pot?_orig_rel = max(tsAccPot?_orig_rel.data);',ic)
  c_eval('acc_pot?_nobg_rel = max(tsAccPot?_nobg_rel.data);',ic)
  
  %% Plot figure of density fred, flux velocity etc
  c_eval('legends_mms{?,1} = sprintf(''mms %g'',?);',ic)
  
  ic = 1;
  nPanels = 7;
  h = irf_plot(nPanels);
  
  if 1 % B
    hca = irf_panel('B');
    set(hca,'ColorOrder',mms_colors('xyz'))
    irf_plot(hca,{gseB1.x,gseB1.y,gseB1.z},'comp')
    hca.YLabel.String = {'B','(nT)'};
    irf_legend(hca,{'x','y','z'},[0.02 0.98])
  end
  if 1 % ne
    hca = irf_panel('ne');
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,ne1,'comp')
    hca.YLabel.String = {'n_e','(cm^{-3})'};    
    hca.YLabel.Interpreter = 'tex';
    irf_legend(hca,{sprintf('n_e^{lb} = %.3f cc',n_lobe)},[0.02 0.98])
    
    %ax2 = axes('Position',get(hca,'Position'));
    %set(ax2,'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required
    %set(ax2,'YAxisLocation','right');
    %set(ax2,'Color','none'); % color of axis

    %yscale = n_lobe;
    %ax2 = gca;
    %ax2.YLim = hca.YLim/yscale;
    %ax2.YLabel.String = 'n_e/n^{lb}';
  end
  if 1 % Ve
    hca = irf_panel('Ve');
    set(hca,'ColorOrder',mms_colors('xyz'))
    irf_plot(hca,{gseVe1.x,gseVe1.y,gseVe1.z},'comp')
    hca.YLabel.String = {'v_e','(km/s)'};
    irf_legend(hca,{'x','y','z'},[0.02 0.98])
    irf_legend(hca,{sprintf('v_e^{lb} = %.0f km/s',cn_eV2v(Tepar_lobe,'eV'))},[0.02 0.1],'k')
  end
  if 1 % Te
    hca = irf_panel('Te');
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,{facTe1.xx,0.5*(facTe1.yy+facTe1.zz)},'comp')
    hca.YLabel.String = {'Te','(eV)'};
    irf_legend(hca,{'||','\perp'},[0.02 0.98])
  end
  if 1 % ne
    hca = irf_panel('ne/n_lb');
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,ne1/n_lobe,'comp')
    hca.YLabel.String = {'n_e/n_e^{lb}',''};    
    hca.YLabel.Interpreter = 'tex';
    irf_legend(hca,{sprintf('n_e^{lb} = %.3f cc',n_lobe)},[0.02 0.98])
  end
  if 1 % Ve
    hca = irf_panel('Ve/vtelbs');
    set(hca,'ColorOrder',mms_colors('xyz'))
    irf_plot(hca,{gseVe1.x/cn_eV2v(Tepar_lobe,'eV'),gseVe1.y/cn_eV2v(Tepar_lobe,'eV'),gseVe1.z/cn_eV2v(Tepar_lobe,'eV')},'comp')
    hca.YLabel.String = {'v_e/v_{te}^{lb}',''};
    irf_legend(hca,{'x','y','z'},[0.02 0.98])
  end
  if 1 % Ve
    hca = irf_panel('ne/n_lb*Ve/vtelbs');
    set(hca,'ColorOrder',mms_colors('xyz'))
    toplot = gseVe1*ne1/cn_eV2v(Tepar_lobe,'eV')/n_lobe;
    irf_plot(hca,{toplot.x,toplot.y,toplot.z},'comp')
    hca.YLabel.String = {'n_ev_e/n_e^{lb}v_{te}^{lb}',''};
    irf_legend(hca,{'x','y','z'},[0.02 0.98])
  end
  c_eval('irf_pl_mark(h(?),tint_lobe,''b'')',1:nPanels)
  c_eval('irf_pl_mark(h(?),tint_sep,''y'')',1:nPanels)
  
  %irf_zoom(h,'x',tint_fred)
  irf_zoom(h,'y')
  h(1).Title.String = sprintf('MMS %g',ic);
  
  %% Plot Liouville mapping of lobe population, based on potential
  iicc = 1;
  load(sprintf('%s/acc_pot_data_event_%g',saveAccPotPath,19),'acc_pot_data')
  c_eval('acc_pot? = acc_pot_data.tsAccPot_nobg_abs{?};,',iicc)
  tt = ef1D1_orig.time;  
  tphi = 8;
  ttrans = 1;
  acc_pot = irf.ts_scalar(tt,1000*(1+1*tanh((tt-tt(1)-tphi)/ttrans)));
  %irf_plot({acc_pot1.abs,acc_pot},'comp') 
  c_eval('[f?_lio,n?_lio] = acc_pot_lio_map(ef1D?_orig.tlim(tint_lobe),acc_pot);',iicc)
  %%
  ic = 1;
  h = irf_plot(5);
  
  if 1 % f red orig.
    hca = irf_panel('f red. orig.');
    c_eval('irf_spectrogram(hca,ef1D?_orig.specrec);',ic)
  end
  if 1 % f red nobg.
    hca = irf_panel('f red. nobg.');
    c_eval('irf_spectrogram(hca,ef1D?_nobg.specrec);',ic)
  end
  if 1 % f red nobg. model
    hca = irf_panel('f red. model');
    c_eval('irf_spectrogram(hca,f?_lio.specrec);',ic)
  end
  linkprop(h(1:3),'CLim')
  
  %cn.print(sprintf('event%g_tsAccPot%g_nobg_abs',event,iic))
  %% Save data to file
  % .mat or .txt (together or separate?)  
  acc_pot_data.event_id = event;
  acc_pot_data.ic = ic;
  acc_pot_data.vix = vix;
  acc_pot_data.bix = bix;
  acc_pot_data.B_lobe = B_lobe;
  acc_pot_data.n_lobe = n_lobe;
  acc_pot_data.n_sheet = n_sheet;
  acc_pot_data.n_sep = n_sep;
  acc_pot_data.n_sep_min = n_sep_min;
  acc_pot_data.be_lobe = be_lobe;
  acc_pot_data.be_sheet = be_sheet;
  acc_pot_data.be_sep = be_sep;
  acc_pot_data.Tepar_lobe = Tepar_lobe;
  acc_pot_data.Tepar_sheet = Tepar_sheet;
  acc_pot_data.Tepar_sep = Tepar_sep;
  acc_pot_data.Teperp_lobe = Teperp_lobe;
  acc_pot_data.Teperp_sheet = Teperp_sheet;
  acc_pot_data.Teperp_sep = Teperp_sep;
  acc_pot_data.bg_nSecondary = nSecondary;
  acc_pot_data.bg_nPhoto = nPhoto;
  acc_pot_data.red_nMC = nMC;
  acc_pot_data.acc_pot_vint = vint;
  acc_pot_data.acc_pot_eint = eint;
  acc_pot_data.acc_pot_resample_timestep = resample_timestep;
  acc_pot_data.acc_pot_minpeakprominence = minpeakprominence;
  acc_pot_data.tsAccPot_orig_abs = tsAccPot_orig_abs;
  acc_pot_data.tsAccPot_nobg_abs = tsAccPot_nobg_abs;
  acc_pot_data.acc_pot_relativeminpeakprominence = relativeminpeakprominence;  
  acc_pot_data.tsAccPot_orig_rel = tsAccPot_orig_rel;
  acc_pot_data.tsAccPot_nobg_rel = tsAccPot_nobg_rel;
  
  save(sprintf('%s/acc_pot_data_event_%g',saveAccPotPath,event),'acc_pot_data')
  %fid = fopen([matlabPath 'acc_pot_statistics.txt'],'a+');
  % event_id
  %save_format = '%f';
  %fid = fclose(fid);
  %fprintf('printed to file %sesw_properties.txt: %s /n',matlabPath,print_str)
end