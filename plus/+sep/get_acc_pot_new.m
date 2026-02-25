% In this new version, use min peak height and min peak prominence

%% Load datastore
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');

doSave = 1;

events = 1:20;
nevents = numel(events);
saveAccPotPath = ['/Users/' localuser '/MATLAB/cn-matlab/+sep/acc_potential_new/'];
printAccPotPath = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/acceleration_potential/'];

%% Run through events
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
  tint_phi_ = tint_phi + 1.5*[-1 1];
  c_eval('eDist?_orig = ePDist?.tlim(tint_phi_);',ic)  
  % Remove background, try first with these default values, mostly for
  % cosmetics anyways
  nSecondary = 5;
  nPhoto = 1;
  c_eval('eDist?_nobg = mms.remove_edist_background(eDist?_orig,''nSecondary'',nSecondary,''Nphotoe_art'',nPhoto);',ic)

  %% Compute reduced electron distribution
  eint = [00 40000];
  lowerelim = 100; % not used when removing background
      
  vgmax = 70000;
  vg = -vgmax:1000:vgmax;
  vg(abs(vg)>70000) = [];

  nMC = 500;
  
  c_eval('tic; ef1D?_orig = eDist?_orig.reduce(''1D'',dmpaB?.resample(eDist?_orig).norm,''scpot'',scPot?.resample(eDist?_orig),''nMC'',nMC,''vg'',vg); toc',ic) % reduced distribution along B
  c_eval('tic; ef1D?_lowe = eDist?_orig.reduce(''1D'',dmpaB?.resample(eDist?_orig).norm,''scpot'',scPot?.resample(eDist?_orig),''lowerelim'',lowerelim,''nMC'',nMC,''vg'',vg); toc',ic) % reduced distribution along B
  c_eval('tic; ef1D?_nobg = eDist?_nobg.reduce(''1D'',dmpaB?.resample(eDist?_nobg).norm,''scpot'',scPot?.resample(eDist?_nobg),''nMC'',nMC,''vg'',vg); toc',ic) % reduced distribution along B

  c_eval('ef1D_orig{?} = ef1D?_orig;',ic)
  c_eval('ef1D_lowe{?} = ef1D?_lowe;',ic)
  c_eval('ef1D_nobg{?} = ef1D?_nobg;',ic)
  
  %% Get acceleration potential
  Te_lobe = (Tepar_lobe + 2*Teperp_lobe)/3;
  vte_lobe = sqrt(2*units.eV*Te_lobe/units.me);
  f_lobe = 1e6*n_lobe/sqrt(pi)/vte_lobe; % sm^-4
  % Decreased interval of tint_phi so that its not always including lobes 
  % or sheet now. Use density and thermal velocity instead.
  
  f_prom = 0.1*f_lobe;
  f_min_peak = 0.2*f_lobe;
  
  eint = [000 40000];
  if any(quadrant([1 4])) % top left or bottom right, inflow is antiparallel to B
    vint = [-100000 -5000];
  else  % top right or bottom left, inflow is parallel to B
    vint = [5000 100000];
  end
  %relativeminpeakprominence = f_prom_abs;
  minpeakheight = f_min_peak;
  minpeakprominence = f_prom;
  resample_timestep = 3;
  % Need to make for-loop for printing, since figure is done within script
  clear tsAccPot_orig tsAccPot_lowe tsAccPot_nobg
  for iic = ic
    % There seems to be no such thing as relativeminpeakprominence
    %% Absolute peak prominence and min peak height
    [tmpTsAccPot,tmppPeakF] = find_acc_pot(ef1D_orig{iic},...
      'eint',eint,...
      'vint',vint,...
      'minpeakheight',minpeakheight,... % fraction of f_lobe
      'minpeakprominence',minpeakprominence,... % smaller fraction of f_lobe
      'resample',resample_timestep);
    % Collect data
    tsAccPot_orig{iic} = tmpTsAccPot;
    peakF_orig{iic} = tmppPeakF;
    % Adjust figure and print
    hcf = gcf;
    hcf.Position = [0 100 700 800];
    h = findobj(hcf.Children,'type','axes');
    c_eval('irf_pl_mark(h(?),tint_phi,''r'');',1:numel(h))
    c_eval('h(?).CLim = log10([0.1*f_lobe f_lobe]);',2:3)
    cn.print(sprintf('event%g_tsAccPot%g_orig',event,iic),'path',printAccPotPath)    
    %%
    %break
    %% Absolute peak prominence and min peak height
    [tmpTsAccPot,tmppPeakF] = find_acc_pot(ef1D_lowe{iic},...
      'eint',eint,...
      'vint',vint,...
      'minpeakheight',minpeakheight,... % fraction of f_lobe
      'minpeakprominence',minpeakprominence,... % smaller fraction of f_lobe
      'resample',resample_timestep);
    % Collect data
    tsAccPot_lowe{iic} = tmpTsAccPot;
    peakF_lowe{iic} = tmppPeakF;        
    % Adjust figure and print
    hcf = gcf;   
    hcf.Position = [0 100 700 800];
    h = findobj(hcf.Children,'type','axes');
    c_eval('h(?).CLim = log10([0.1*f_lobe f_lobe]);',2:3)
    cn.print(sprintf('event%g_tsAccPot%g_lowe',event,iic),'path',printAccPotPath)
    
    %% Absolute peak prominence and min peak height
    [tmpTsAccPot,tmppPeakF] = find_acc_pot(ef1D_nobg{iic},...
      'eint',eint,...
      'vint',vint,...
      'minpeakheight',minpeakheight,... % fraction of f_lobe
      'minpeakprominence',minpeakprominence,... % smaller fraction of f_lobe
      'resample',resample_timestep);
    % Collect data
    tsAccPot_nobg{iic} = tmpTsAccPot;
    peakF_nobg{iic} = tmppPeakF;        
    % Adjust figure and print
    hcf = gcf;   
    hcf.Position = [0 100 700 800];
    h = findobj(hcf.Children,'type','axes');
    c_eval('h(?).CLim = log10([0.1*f_lobe f_lobe]);',2:3)
    cn.print(sprintf('event%g_tsAccPot%g_nobg',event,iic),'path',printAccPotPath)
    
  end
  
  % Already doing this above now
%   c_eval('tsAccPot_orig{?} = tsAccPot?_orig;',ic) 
%   c_eval('tsAccPot_lowe{?} = tsAccPot?_lowe;',ic)
%   c_eval('tsAccPot_nobg{?} = tsAccPot?_nobg;',ic)  
%           
%   % From timeseries, pick out maximal value
%   % Can do this when reloading data
%   c_eval('acc_pot?_orig = max(tsAccPot?_orig.data);',ic)
%   c_eval('acc_pot?_lowe = max(tsAccPot?_lowe.data);',ic)
%   c_eval('acc_pot?_nobg = max(tsAccPot?_nobg.data);',ic)  
  
  %% Plot acceleration potentials for all spacecraft
  c_eval('legends_mms{?,1} = sprintf(''mms %g'',?);',ic)
  
  h = irf_plot(3);
  
  if 1 % orig
    hca = irf_panel('acc. pot. orig.');
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,tsAccPot_orig,'comp')
    hca.YLabel.String = {'Acc. pot. orig.','(eV)'};
    irf_legend(hca,legends_mms,[0.02 0.98])
  end
  if 1 % lowe
    hca = irf_panel('acc. pot. orig. rel.');
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,tsAccPot_lowe,'comp')
    hca.YLabel.String = {'Acc. pot. lowe.','(eV)'};
    irf_legend(hca,legends_mms,[0.02 0.98])
  end
  if 1 % nobg
    hca = irf_panel('acc. pot. nobg.');
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_plot(hca,tsAccPot_nobg,'comp')
    hca.YLabel.String = {'Acc. pot. nobg.','(eV)'};
    irf_legend(hca,legends_mms,[0.02 0.98])
  end
  irf_zoom(h,'x',tint_phi_)
  
  %% Plot 1D phase space density that gave maximum beam energy
  if numel(ic) == 3
    colors = mms_colors('123');
  elseif  numel(ic) == 4
    colors = mms_colors('1234');
  end
  h = setup_subplots(3,1);
  isub = 1;
  if 1 % orig
    hca = h(isub); isub = isub + 1;
    ff = peakF_orig;
    plot(hca,ff{1}.v*1e-3,ff{1}.v*0+f_min_peak)
    hold(hca,'on')
    set(hca,'ColorOrder',colors)
    for iic = ic      
      plot(hca,ff{iic}.v*1e-3,ff{iic}.f_orig/f_lobe,'color',colors(iic,:))%,peakF_nobg{iic}.vbeam,peakF_nobg{iic}.fbeam,'o')
      errorbar(hca,ff{iic}.vbeam*1e-3,ff{iic}.fbeam/f_lobe,f_prom/f_lobe,0,'color',colors(iic,:))
    end    
    hold(hca,'off')
    hca.XLim = ff{1}.v([1 end])*1e-3;
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'f/f_{lobe}';
    irf_legend(hca,{'f_{orig}';sprintf('f_{lobe} = %g s/m^4',f_lobe)},[0.02 0.98])
  end
  if 1 % lowe
    hca = h(isub); isub = isub + 1;
    ff = peakF_lowe;
    plot(hca,ff{1}.v*1e-3,ff{1}.v*0+f_min_peak)
    hold(hca,'on')
    set(hca,'ColorOrder',colors)
    for iic = ic      
      plot(hca,ff{iic}.v*1e-3,ff{iic}.f_orig/f_lobe,'color',colors(iic,:))%,peakF_nobg{iic}.vbeam,peakF_nobg{iic}.fbeam,'o')
      errorbar(hca,ff{iic}.vbeam*1e-3,ff{iic}.fbeam/f_lobe,f_prom/f_lobe,0,'color',colors(iic,:))
    end    
    hold(hca,'off')
    hca.XLim = ff{1}.v([1 end])*1e-3;
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'f/f_{lobe}';
    irf_legend(hca,{'f_{lowe}';sprintf('f_{lobe} = %g s/m^4',f_lobe)},[0.02 0.98])
  end
  if 1 % nobg
    hca = h(isub); isub = isub + 1;
    ff = peakF_nobg;
    plot(hca,ff{1}.v*1e-3,ff{1}.v*0+f_min_peak)
    hold(hca,'on')
    set(hca,'ColorOrder',colors)
    for iic = ic      
      plot(hca,ff{iic}.v*1e-3,ff{iic}.f_orig/f_lobe,'color',colors(iic,:))%,peakF_nobg{iic}.vbeam,peakF_nobg{iic}.fbeam,'o')
      errorbar(hca,ff{iic}.vbeam*1e-3,ff{iic}.fbeam/f_lobe,f_prom/f_lobe,0,'color',colors(iic,:))
    end    
    hold(hca,'off')
    hca.XLim = ff{1}.v([1 end])*1e-3;
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'f/f_{lobe}';
    irf_legend(hca,{'f_{nobg}';sprintf('f_{lobe} = %g s/m^4',f_lobe)},[0.02 0.98])
  end
  ylimmax = 0;
  for ip = 1:3
    h(ip).XGrid = 'on';
    h(ip).YGrid = 'on';
    ylimmax = max([ylimmax h(ip).YLim(2)]);    
  end
  hlinks = linkprop(h,{'YLim'});
  hlinks.Targets(1).YLim(2) = ylimmax;
  compact_panels(0.01)
  
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
  acc_pot_data.acc_pot_minpeakheight = minpeakheight;
  acc_pot_data.tsAccPot_orig = tsAccPot_orig;
  acc_pot_data.tsAccPot_lowe = tsAccPot_lowe;
  acc_pot_data.tsAccPot_nobg = tsAccPot_nobg;
  
  acc_pot_data.ef1D_orig = ef1D_orig;
  acc_pot_data.ef1D_lowe = ef1D_lowe;
  acc_pot_data.ef1D_nobg = ef1D_nobg;
  
  acc_pot_data.peakF_orig = peakF_orig;
  acc_pot_data.peakF_lowe = peakF_lowe;
  acc_pot_data.peakF_nobg = peakF_nobg;
  
  acc_pot_data.fe_lobe = f_lobe;
  acc_pot_data.f_prom = f_prom;
  acc_pot_data.f_min_peak = f_min_peak;
  
  
  if doSave
    save(sprintf('%s/acc_pot_data_event_%g',saveAccPotPath,event),'acc_pot_data')
  end
  %fid = fopen([matlabPath 'acc_pot_statistics.txt'],'a+');
  % event_id
  %save_format = '%f';
  %fid = fclose(fid);
  %fprintf('printed to file %sesw_properties.txt: %s /n',matlabPath,print_str)
end