% Rerun of get_acc_pot_new to test different peak heights etc. Reduced
% dists are saved to .mat file.

%% Load datastore
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');

units = irf_units;
doSave = 1;
doPrint = 1;

events = 1:20;
nevents = numel(events);
loadAccPotPath = ['/Users/' localuser '/MATLAB/cn-matlab/+sep/acc_potential_new/'];
saveAccPotPath = ['/Users/' localuser '/MATLAB/cn-matlab/+sep/acc_potential_new_rerun/'];
printAccPotPath = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/acceleration_potential/'];

allEvents = 1:nevents;
okEvents = [1 5 8];
badEvents = setdiff(allEvents,okEvents);

%% Run through events and check whats in there
for ievent = badEvents
  event = events(ievent); % set event id
  sep.get_tints; % get tints for that event id
 
  ic = 1:4;
    
  %% Load saved acc data
  load(sprintf('%s/acc_pot_data_event_%g',saveAccPotPath,event),'acc_pot_data'); % acc_pot_data
  tint_dist = acc_pot_data.ef1D_orig{1}.time([1 end]);
  fields = fieldnames(acc_pot_data);
  for ifield = 1:numel(fields)
    eval(sprintf('%s = acc_pot_data.(''%s'');',fields{ifield},fields{ifield}))
  end
  f_lobe = fe_lobe;

  
  
  if 0 % Plot acceleration potentials for all spacecraft
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
    irf_zoom(h,'x',tint_dist)
  end
  if 1 % Plot 1D phase space density that gave maximum beam energy
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
      ff = acc_pot_data.peakF_nobg;
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
    %hlinks.Targets(1).YLim(2) = ylimmax;
    hlinks.Targets(1).YLim = [0 0.5];
    compact_panels(0.01)
    if doPrint
      cn.print(sprintf('iev=%02.0f_1D_dists_ylim_0_05',ievent),'path',printAccPotPath)
    end
  end
    
  if 0 % re-collect data
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
  end
  %pause
end

%% Run through events and redo analysis
climfraction = [0.01 1];
for ievent = badEvents(7:end)%1:nevents
  event = events(ievent); % set event id
  sep.get_tints; % get tints for that event id
  tint_phi_ = tint_phi + 1.5*[-1 1];
  %% Load saved acc data
  load(sprintf('%s/acc_pot_data_event_%g',loadAccPotPath,event),'acc_pot_data'); % acc_pot_data
  tint_dist = acc_pot_data.ef1D_orig{1}.time([1 end]);
  fields = fieldnames(acc_pot_data);
  for ifield = 1:numel(fields)
    eval(sprintf('%s = acc_pot_data.(''%s'');',fields{ifield},fields{ifield}))
  end
  f_lobe = fe_lobe;
  
  %% Pick out reconnection quadrant: tail, Earth to the left
  quadrant = ones(2,2);  
  if vix > 0; quadrant(:,2) = 0; else quadrant(:,1) = 0; end      
  if bix > 0; quadrant(2,:) = 0; else quadrant(1,:) = 0; end
  
  %% Get acceleration potential
  Te_lobe = (Tepar_lobe + 2*Teperp_lobe)/3;
  vte_lobe = sqrt(2*units.eV*Te_lobe/units.me);
  f_lobe = 1e6*n_lobe/sqrt(pi)/vte_lobe; % sm^-4
  % Decreased interval of tint_phi so that its not always including lobes t
  % or sheet now. Use density and thermal velocity instead.
  
  clim = climfraction*f_lobe;
  f_prom = 0.05*f_lobe;
  f_min_peak = 0.1*f_lobe;
  
  eint = [000 40000];
  if any(quadrant([1 4])) % top left or bottom right, inflow is antiparallel to B
    vint = [-100000 -3000];
  else  % top right or bottom left, inflow is parallel to B
    vint = [3000 100000];
  end
  %relativeminpeakprominence = f_prom_abs;
  minpeakheight = f_min_peak;
  minpeakprominence = f_prom;
  resample_timestep = 3;
  % Need to make for-loop for printing, since figure is done within script
  clear tsAccPot_orig tsAccPot_lowe tsAccPot_nobg ...
    peakF_orig peakF_lowe peakF_nobg
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
    c_eval('h(?).CLim = log10(clim);',2:3)
    cn.print(sprintf('event%g_tsAccPot%g_orig',event,iic),'path',printAccPotPath)    

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
    c_eval('h(?).CLim = log10(clim);',2:3)
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
    c_eval('h(?).CLim = log10(clim);',2:3)
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
    %irf_legend(hca,subsref(ff.time.utc,struct('type','()','subs',{{1:23}})),[0.98 0.98])
    times = cellfun(@(x) subsref(x.time.utc,struct('type','()','subs',{{1:23}})),ff,'UniformOutput',0);    
    irf_legend(hca,times',[0.98 0.98])    
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
    times = cellfun(@(x) subsref(x.time.utc,struct('type','()','subs',{{1:23}})),ff,'UniformOutput',0);    
    irf_legend(hca,times',[0.98 0.98])
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
    times = cellfun(@(x) subsref(x.time.utc,struct('type','()','subs',{{1:23}})),ff,'UniformOutput',0);    
    irf_legend(hca,times',[0.98 0.98])
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
  acc_pot_data.bg_nSecondary = bg_nSecondary;
  acc_pot_data.bg_nPhoto = bg_nPhoto;
  acc_pot_data.red_nMC = red_nMC;
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