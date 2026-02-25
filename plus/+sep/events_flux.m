%% Load datastore
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');

events = 1:20;
nevents = numel(events);
saveAccPotPath = ['/Users/' localuser '/MATLAB/cn-matlab/+sep/acc_potential/'];
printAccPotPath =  ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/flux/'];

for ievent = 11:nevents
  event = events(ievent); % set event id
  sep.get_tints; % get tints for that event id
  
  %% Load data  
  ic = 1:4;
  ic_pdist = 1;
  units = irf_units;
  % Core data
  c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
  if isempty(dmpaB4); ic = 1:3; end
  c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
  c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
  c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic_pdist)  
    
  % Auxillary data, for context
  c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
  c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
  c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
  c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
  c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);

  % Auxillary quantities
  c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
  c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)  
  c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
  c_eval('beta?e = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
  c_eval('beta?i = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
  c_eval('beta? = beta?i + beta?e.resample(beta?i);',ic)
  c_eval('gseFlux? = ne?*gseVe?;',ic);
  c_eval('gseFlux?par = ne?*gseVe?par;',ic);
  c_eval('vte?par = (2*units.eV*facTe?.xx/units.me).^0.5;',ic)
  c_eval('vte?perp = (2*units.eV*(1/2)*(facTe?.yy+facTe?.zz)/units.me).^0.5;',ic)
  c_eval('vte? = (2*units.eV*(1/3)*(facTe?.xx+facTe?.yy+facTe?.zz)/units.me).^0.5;',ic)
  
  units = irf_units;
  % Core data   
  c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
  cad_resamp = 5;
  timeline = ne1.time(1:cad_resamp:end);
  c_eval('ne?_resamp = ne?.resample(timeline);',ic);
  
  %% Pick out parameters in the lobe, sheet, and separatrix intervals
  ic_orig = ic;
  ic = 1; 

  c_eval('n_lobe(?) = mean(ne?.tlim(tint_lobe).data,1);',ic)
  c_eval('n_sheet(?) = mean(ne?.tlim(tint_sheet).data,1);',ic)
  c_eval('n_sep(?) = mean(ne?.tlim(tint_sep).data,1);',ic)
  c_eval('n_sep_min(?) = min(ne?.tlim(tint_sep).data);',ic)
%  c_eval('n_sep_min(?) = min(ne?.resample(.tlim(tint_sep).data);',ic)

  colors = mms_colors('xzy');
  info_color = [0 0 0; colors; colors; colors; 0 0 0; 0 0 0; 0 0 0; colors];

  tint_day_utc = tint.utc; tint_day_utc = tint_day_utc(1,1:10);
  tint_lobe_utc = tint_lobe.utc; tint_lobe_utc = tint_lobe_utc(:,12:23);
  tint_sheet_utc = tint_sheet.utc; tint_sheet_utc = tint_sheet_utc(:,12:23);
  tint_sep_utc = tint_sep.utc; tint_sep_utc = tint_sep_utc(:,12:23);
  tint_phi_utc = tint_phi.utc; tint_phi_utc = tint_phi_utc(:,12:23);
  
  ic = ic_orig;     
  
  %% Pick out parameters in the lobe, sheet, and separatrix intervals
  ic_orig = ic;  
  c_eval('B_lobe(?) = mean(gseB?.abs.tlim(tint_lobe).data,1);',ic)

  c_eval('n_lobe(?) = mean(ne?.tlim(tint_lobe).data,1);',ic)
  c_eval('n_sheet(?) = mean(ne?.tlim(tint_sheet).data,1);',ic)
  c_eval('n_sep(?) = mean(ne?.tlim(tint_sep).data,1);',ic)
  c_eval('n_sep_min(?) = min(ne?.tlim(tint_sep).data);',ic)
  %c_eval('n_sep_data(?) = ne?.tlim(tint_sep);',ic)

  c_eval('be_lobe(?) = mean(beta?e.tlim(tint_lobe).data,1);',ic)
  c_eval('be_sheet(?) = mean(beta?e.tlim(tint_sheet).data,1);',ic)
  c_eval('be_sep(?) = mean(beta?e.tlim(tint_sep).data,1);',ic)

  c_eval('Tepar_lobe(?) = mean(facTe?.tlim(tint_lobe).data(:,1,1),1);',ic)
  c_eval('Tepar_sheet(?) = mean(facTe?.tlim(tint_sheet).data(:,1,1),1);',ic)
  c_eval('Tepar_sep(?) = mean(facTe?.tlim(tint_sep).data(:,1,1),1);',ic)

  c_eval('Teperp_lobe(?)  = 0.5*mean(facTe?.tlim(tint_lobe).data(:,2,2)+facTe?.tlim(tint_lobe).data(:,3,3),1);',ic)
  c_eval('Teperp_sheet(?) = 0.5*mean(facTe?.tlim(tint_sheet).data(:,2,2)+facTe?.tlim(tint_sheet).data(:,3,3),1);',ic)
  c_eval('Teperp_sep(?)   = 0.5*mean(facTe?.tlim(tint_sep).data(:,2,2)+facTe?.tlim(tint_sep).data(:,3,3),1);',ic)

  c_eval('vteperp_lobe(?) = sqrt(2*units.eV*Teperp_lobe(?)/units.me)*1e-3;',ic)
  c_eval('vtepar_lobe(?) = sqrt(2*units.eV*Tepar_lobe(?)/units.me)*1e-3;',ic)
  c_eval('vte_lobe(?) = sqrt(2*units.eV*(1/3)*(2*Teperp_lobe(?)+Tepar_lobe(?))/units.me)*1e-3;',ic)
  
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

  %% Compute reduced electron distribution
  % %% Remove photoelectron and secondary electron background, for tint_phi
  %tint_fred = tint_phi;%tint_phi;
  tint_phi_ = tint_fred;
  tint_phi_ = tint_flow + [-2 2];
  c_eval('eDist?_orig = ePDist?.tlim(tint_phi_);',ic_pdist)  
  % Remove background, try first with these default values, mostly for
  % cosmetics anyways
  %nSecondary = 5;
  %nPhoto = 1;
  %c_eval('eDist?_nobg = mms.remove_edist_background(eDist?_orig,''nSecondary'',nSecondary,''Nphotoe_art'',nPhoto);',ic)

  
  eint = [00 40000];
  lowerelim = 000; % not used when removing background
      
  vgmax = 70000;
  vg = -vgmax:1000:vgmax;
  vg(abs(vg)>70000) = [];

  nMC = 200;
  
  c_eval('tic; ef1D?_orig = eDist?_orig.reduce(''1D'',dmpaB?.resample(eDist?_orig).norm,''scpot'',scPot?.resample(eDist?_orig),''lowerelim'',lowerelim,''nMC'',nMC,''vg'',vg); toc',ic_pdist) % reduced distribution along B
%  c_eval('tic; ef1D?_nobg = eDist?_nobg.reduce(''1D'',dmpaB?.resample(eDist?_nobg).norm,''scpot'',scPot?.resample(eDist?_nobg),''lowerelim'',lowerelim,''nMC'',nMC,''vg'',vg); toc',ic) % reduced distribution along B

  % Liouville mapping  
  % get map from two_dim_dependence_of_liouville_density...
  c_eval('[ts_vtr?,ts_vpsi?] = intersection_2d(sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),Nmap,NVmap,ne?/n_lobe(?),gseFlux?par.abs/(n_lobe(?)*vte_lobe(?)),''plot'',0,''mesh'',0);',ic)
  c_eval('[ts_vtr?sep,ts_vpsi?sep] = intersection_2d(sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),Nmap,NVmap,ne?.tlim(tint_sep)/n_lobe(?),gseFlux?par.tlim(tint_sep).abs/(n_lobe(?)*vte_lobe(?)),''plot'',0,''mesh'',0);',ic)
  
  c_eval('[ts_vtr?flow,ts_vpsi?flow] = intersection_2d(sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),Nmap,NVmap,ne?.tlim(tint_flow)/n_lobe(?),gseFlux?par.tlim(tint_flow).abs/(n_lobe(?)*vte_lobe(?)),''plot'',0,''mesh'',0);',ic)
  
  % multiply flux with 2
  c_eval('[ts_vtr?sep_2,ts_vpsi?sep_2] =   intersection_2d(sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),Nmap,NVmap,ne?.tlim(tint_sep)/n_lobe(?),2*gseFlux?par.tlim(tint_sep).abs/(n_lobe(?)*vte_lobe(?)),''plot'',0,''mesh'',0);',ic)
  c_eval('[ts_vtr?flow_2,ts_vpsi?flow_2] = intersection_2d(sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),Nmap,NVmap,ne?.tlim(tint_flow)/n_lobe(?),2*gseFlux?par.tlim(tint_flow).abs/(n_lobe(?)*vte_lobe(?)),''plot'',0,''mesh'',0);',ic)
   
  %ts_fred1flow = get_lioumap_dist(n_lobe(1),te_lobe(1),ts_vtr1flow_2,ts_vpsi1flow_2);
  ts_fred1flow = get_lioumap_dist(n_lobe(1),Tepar_lobe(1),ts_vtr1flow,ts_vpsi1flow);  
  ts_fred1flow_2 = get_lioumap_dist(n_lobe(1),Tepar_lobe(1),ts_vtr1flow_2,ts_vpsi1flow_2);  
  if 1
    %%
    if 0
    ic_=ic;
    ic = 1;
    h = irf_plot(4);
    
    if 1 % vtr
      hca = irf_panel('vtr');
      set(hca,'ColorOrder',mms_colors('1234'))
      irf_plot(hca,ts_vtr1flow,'comp');
      hca.YLabel.String = {'(\psi/T_0)^{1/2}',''};
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end
    if 1 % vpsi
      hca = irf_panel('vpsi');
      set(hca,'ColorOrder',mms_colors('1234'))
      irf_plot(hca,ts_vpsi1flow,'comp');
      hca.YLabel.String = {'v_\psi/v_{te}^{lb}',''};
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end       
    if 1 % fred
      hca = irf_panel('reduced fe');
      set(hca,'ColorOrder',mms_colors('1234'))
      c_eval('specrec = ef1D?_orig.specrec(''10^3 km/s'');',ic)  
      flim = 10^-7;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end    
    if 1 % fred
      hca = irf_panel('reduced fe mapped');
      set(hca,'ColorOrder',mms_colors('1234'))
      c_eval('specrec = ts_fred?flow.specrec(''10^3 km/s'');',ic)  
      flim = 10^-10;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(...)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
      hca.YLim = get(irf_panel('reduced fe'),'ylim');
      hca.CLim = get(irf_panel('reduced fe'),'clim');
    end
    irf_plot_axis_align(h)
    irf_zoom(h,'x',ts_vpsi1.time([1 end]))
    ic = ic_;
    end
    %%
    ic_=ic;
    ic = 1;
    printdir = [printAccPotPath 'flux_therm/event' num2str(ievent) '/'];
    mkdir(printdir);

    c_eval('nlev = ne?.tlim(tint_flow)/n_lobe(?);',ic)
    c_eval('nvlev = gseFlux?par.tlim(tint_flow)/(n_lobe(?)*vte_lobe(?));',ic)
      
    for it = 1:ts_fred1flow.length
      hca = subplot(1,2,1);
      if nvlev(it).data > 0
        plot(hca,...
          ef1D1_orig.tlim(tint_flow).depend{1}(1,:)*1e-3,ef1D1_orig.tlim(tint_flow).data(it,:),...          
          -ts_fred1flow.depend{1}(1,:)*1e-3,ts_fred1flow.data(it,:)...,...
          ...ts_fred1flow_2.depend{1}(1,:)*1e-3,ts_fred1flow_2.data(it,:)
          ); 
      else        
        plot(hca,...
          ef1D1_orig.tlim(tint_flow).depend{1}(1,:)*1e-3,ef1D1_orig.tlim(tint_flow).data(it,:),...
          ts_fred1flow.depend{1}(1,:)*1e-3,ts_fred1flow.data(it,:)...          
          ...ts_fred1flow_2.depend{1}(1,:)*1e-3,ts_fred1flow_2.data(it,:)
          ); 
      end
      
      title(hca,sprintf('it=%g',it)); 
      hca.XLim = 50e3*[-1 1]*1e-3;
      hca.XTick = [-50:10:50];
      hca.YLim = [0 max(ts_fred1flow.data(:))*1.4];
      hca.XLabel.String = 'v_{e,||} (10^3 km/s)';
      hca.YLabel.String = 'f_{e}(v_{||}) (s/m^4)';
      hca.XGrid = 'on';
      hca.YGrid = 'on';
      
      hca = subplot(1,2,2);
      colors = pic_colors('matlab');
      [h__,hcc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),Nmap,0.1:0.1:1,'k');
      clabel(h__,hcc);
      hold(hca,'on')
      [h__,hcc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),NVmap,0.2:0.2:3,'r');    
      clabel(h__,hcc);    
      irf_legend(hca,{sprintf('n = %.3f',ne1.tlim(tint_flow).data(it))},[0.02 0.98])
            
      plot(hca,ts_vtr1flow(it).data,ts_vpsi1flow(it).data,'o','color',colors(2,:))
      
      % plot lines 
      [h__,hcc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),Nmap,nlev(it).data*[1 1],'k');    
      hcc.LineWidth = 1.5;
      [h__,hcc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),NVmap,nvlev(it).abs.data*[1 1],'r');    
      hcc.LineWidth = 1.5;
      
      %plot(hca,ts_vtr1flow_2(it).data,ts_vpsi1flow_2(it).data,'o','color',colors(3,:))
      hold(hca,'off')
      hca.XLabel.String = '(\psi/T_{e}^{lb})^{1/2}';
      hca.YLabel.String = 'v_{\psi}/v_{te}^{lb}';
      
      title(hca,sprintf('time = %s',ts_fred1flow(it).time.utc)); 
    
      ic = ic_;
      
      cn.print(sprintf('flux_mms1_event%g_it%g',event,it),'path',printdir)
    end
  end
  
  
  %%
  % Make figures
  if 0 % singel sc, mapping
    %% 
    ic = 1;
    nPanels = 7;
    nRows = 2;
    nCols = 1;
    [h,h2] = initialize_combined_plot(nPanels,nRows,nCols,0.6,'vertical'); 
    h_add = h;
    h_add(:) = [];
    if 1 % fred
      hca = irf_panel('reduced fe');
      set(hca,'ColorOrder',mms_colors('1234'))
      c_eval('specrec = ef1D?_orig.specrec;',ic)  
      flim = 10^-7;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end
    if 1 % B
      hca = irf_panel('B');
      set(hca,'ColorOrder',mms_colors('xyz'))
      c_eval('irf_plot(hca,gseB?);',ic)
      hca.YLabel.String = {'B','(nT)'};
      irf_legend(hca,{'x','y','z'},[0.02 0.98])
    end
    if 1 % n
      hca = irf_panel('n');
      set(hca,'ColorOrder',mms_colors('xyz'))
      c_eval('irf_plot(hca,ne?);',ic)
      hca.YLabel.String = {'n','(cm^{-3})'};
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,{'x','y','z'},[0.02 0.98])

      c_eval('yscale = n_lobe(?);',ic)
      ax2 = axes('Position',get(hca,'Position'));
      set(ax2,'XAxisLocation','bottom','xtick',[]); % remove 'xtick' if xticks required
      set(ax2,'YAxisLocation','right');
      set(ax2,'Color','none'); % color of axis
      set(ax2,'XColor','k','YColor','r'); % color of axis lines and numbers
      ax2.YLim = hca.YLim/yscale;
      ax2.YGrid = 'on';
      ax2.YLabel.String = {'n/n^{lb}',''};
      h_add(end+1) = ax2;

      irf_legend(hca,{['n^{lb} = ' sprintf('%.3f cc',yscale)]},[0.02 0.98],'k')
    end
    if 0 % v xyz
      hca = irf_panel('v');
      set(hca,'ColorOrder',mms_colors('xyz'))
      c_eval('irf_plot(hca,gseVe?);',ic)
      hca.YLabel.String = {'v_e','(km/s)'};
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,{'x','y','z'},[0.02 0.98])

      c_eval('yscale = vte_lobe(?);',ic)
      ax2 = axes('Position',get(hca,'Position'));
      set(ax2,'XAxisLocation','bottom','xtick',[]); % remove 'xtick' if xticks required
      set(ax2,'YAxisLocation','right');
      set(ax2,'Color','none'); % color of axis
      set(ax2,'XColor','k','YColor','r'); % color of axis lines and numbers
      ax2.YLim = hca.YLim/yscale;    
      ax2.YGrid = 'on';
      ax2.YLabel.String = {'v_e/v_{te}^{lb}',''};
      h_add(end+1) = ax2;

      irf_legend(hca,{['v_{te}^{lb} = ' sprintf('%.0f km/s',yscale)]},[0.02 0.98],'k')
    end
    if 1 % v par
      hca = irf_panel('v');
      set(hca,'ColorOrder',mms_colors('xyz'))
      c_eval('irf_plot(hca,gseVe?par);',ic)
      hca.YLabel.String = {'v_{e,||}','(km/s)'};
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,{'x','y','z'},[0.02 0.98])

      c_eval('yscale = vte_lobe(?);',ic)
      ax2 = axes('Position',get(hca,'Position'));
      set(ax2,'XAxisLocation','bottom','xtick',[]); % remove 'xtick' if xticks required
      set(ax2,'YAxisLocation','right');
      set(ax2,'Color','none'); % color of axis
      set(ax2,'XColor','k','YColor','r'); % color of axis lines and numbers
      ax2.YLim = hca.YLim/yscale;    
      ax2.YGrid = 'on';
      ax2.YLabel.String = {'v_e/v_{te}^{lb}',''};
      h_add(end+1) = ax2;

      irf_legend(hca,{['v_{te}^{lb} = ' sprintf('%.0f km/s',yscale)]},[0.02 0.98],'k')
    end
    if 0 % nv xyz
      hca = irf_panel('nv');
      set(hca,'ColorOrder',mms_colors('xyz'))
      c_eval('irf_plot(hca,gseFlux?);',ic)
      hca.YLabel.String = {'nv','(kms^{-1}cm^{-3})'};
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,{'x','y','z'},[0.02 0.98])

      c_eval('yscale = n_lobe(?)*vte_lobe(?);',ic)
      ax2 = axes('Position',get(hca,'Position'));
      set(ax2,'XAxisLocation','bottom','xtick',[]); % remove 'xtick' if xticks required
      set(ax2,'YAxisLocation','right');
      set(ax2,'Color','none'); % color of axis
      set(ax2,'XColor','k','YColor','r'); % color of axis lines and numbers
      ax2.YLim = hca.YLim/yscale;
      ax2.YGrid = 'on';
      ax2.YLabel.String = {'nv/n^{lb}v_{te}^{lb}',''};
      h_add(end+1) = ax2;
    end
    if 1 % nv par
      hca = irf_panel('nv');
      set(hca,'ColorOrder',mms_colors('xyz'))
      c_eval('irf_plot(hca,gseFlux?par);',ic)
      hca.YLabel.String = {'nv_{||}','(kms^{-1}cm^{-3})'};
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,{'x','y','z'},[0.02 0.98])

      c_eval('yscale = n_lobe(?)*vte_lobe(?);',ic)
      ax2 = axes('Position',get(hca,'Position'));
      set(ax2,'XAxisLocation','bottom','xtick',[]); % remove 'xtick' if xticks required
      set(ax2,'YAxisLocation','right');
      set(ax2,'Color','none'); % color of axis
      set(ax2,'XColor','k','YColor','r'); % color of axis lines and numbers
      ax2.YLim = hca.YLim/yscale;
      ax2.YGrid = 'on';
      ax2.YLabel.String = {'nv/n^{lb}v_{te}^{lb}',''};
      h_add(end+1) = ax2;
    end
    if 0 % n/nlb
      hca = irf_panel('n/nlb');
      set(hca,'ColorOrder',mms_colors('xyz'))
      c_eval('irf_plot(hca,ne?/n_lobe(?));',ic)
      hca.YLabel.String = {'n/n^{lb}',''};
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,{'x','y','z'},[0.02 0.98])
    end
    if 0 % v
      hca = irf_panel('v/vtelb');
      set(hca,'ColorOrder',mms_colors('xyz'))
      c_eval('irf_plot(hca,gseVe?/vte_lobe(?));',ic)
      hca.YLabel.String = {'v_e/v_{te}^{lb}',''};
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,{'x','y','z'},[0.02 0.98])
    end
    if 0 % nv
      hca = irf_panel('nv/nlbvtelb');
      set(hca,'ColorOrder',mms_colors('xyz'))
      c_eval('irf_plot(hca,gseFlux?/n_lobe(?)/vte_lobe(?));',ic)
      hca.YLabel.String = {'nv/n^{lb}v_{te}^{lb}',''};
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,{'x','y','z'},[0.02 0.98])
    end

    if 1 % ..
      hca = irf_panel('..');
      set(hca,'ColorOrder',mms_colors('xyz'))
      c_eval('irf_plot(hca,ts_vtr?flow);',ic)
      hca.YLabel.String = {'(\psi/T_{e}^{lb})^{1/2}',''};    
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,{'x','y','z'},[0.02 0.98])
    end
    if 1 % ...
      hca = irf_panel('...');
      set(hca,'ColorOrder',mms_colors('xyz'))
      c_eval('irf_plot(hca,ts_vpsi?flow);',ic)
      hca.YLabel.String = {'v_{\psi}/v_{te}^{lb}',''};
      hca.YLabel.Interpreter = 'tex';
      %irf_legend(hca,{'x','y','z'},[0.02 0.98])
    end

    c_eval('irf_pl_mark(h(?),tint_sep,''y'')',2:nPanels)
    c_eval('irf_pl_mark(h(?),tint_flow,''g'')',2:nPanels)
    c_eval('irf_pl_mark(h(?),tint_sheet,''r'')',2:nPanels)
    c_eval('irf_pl_mark(h(?),tint_lobe,''b'')',2:nPanels)
    irf_zoom(h,'x',tint_fred)
    irf_plot_axis_align([h h_add])

    for ii = 1:numel(h_add)    
      h_add(ii).YLabel.Position(1) = 1.1;
    end


    hca = h2(1);
    [h__,hcc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),Nmap,0.1:0.1:1,'k');
    clabel(h__,hcc);
    hold(hca,'on')
    [h__,hcc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),NVmap,0.2:0.2:3,'r');    
    clabel(h__,hcc);
    plot(hca,ts_vtr1.data,ts_vpsi1.data,'ko')
    plot(hca,ts_vtr1flow.data,ts_vpsi1flow.data,'g*')
    plot(hca,ts_vtr1sep.data,ts_vpsi1sep.data,'y*')
    hold(hca,'off')
    hca.XLabel.String = '(\psi/T_{e}^{lb})^{1/2}';
    hca.YLabel.String = 'v_{\psi}/v_{te}^{lb}';

  
    hca = h2(2);
    scatx = ne1.data/n_lobe(1);
    scaty = gseFlux1par.abs.data/(n_lobe(1)*vte_lobe(1));
    hscat = scatter(hca,scatx,scaty,'.');
    hold(hca,'on')
    scatx = ne1.tlim(tint_flow).data/n_lobe(1);
    scaty = gseFlux1par.tlim(tint_flow).abs.data/(n_lobe(1)*vte_lobe(1));
    hscat = scatter(hca,scatx,scaty,'ok');
    hscat = scatter(hca,scatx,scaty,'.g');
    scatx = ne1.tlim(tint_sep).data/n_lobe(1);
    scaty = gseFlux1par.tlim(tint_sep).abs.data/(n_lobe(1)*vte_lobe(1));
    hscat = scatter(hca,scatx,scaty,'ok');
    hscat = scatter(hca,scatx,scaty,'.y');
    hold(hca,'off')
    hca.XLabel.String = 'n_e/n^{lb}';
    hca.YLabel.String = 'n_ev_{e||}/n^{lb}v_{te}^{lb}';
    hca.Box = 'on';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    %cn.print('flux_mms1','path',printAccPotPath)
    cn.print(sprintf('tint%s-%s_flux_mms1_map_event%g',tint_sep_utc(1,:),tint_sep_utc(2,:),events(ievent)),'path',[printAccPotPath 'flux_map/'])
  end
  if 0 % all three/four together  
    nPanels = 5;
    h = irf_plot(nPanels); 
    h_add = h;
    h_add(:) = [];
    if 1 % fred mms1
      hca = irf_panel('reduced fe mms1');
      set(hca,'ColorOrder',mms_colors('1234'))
      specrec = ef1D1_orig.specrec('10^3 km/s'); 
      flim = 10^-7;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end
    if 0 % fred mms2
      hca = irf_panel('reduced fe mms2');
      set(hca,'ColorOrder',mms_colors('1234'))
      specrec = ef1D2_orig.specrec('10^3 km/s'); 
      flim = 10^-7;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end
    if 0 % fred mms1
      hca = irf_panel('reduced fe mms3');
      set(hca,'ColorOrder',mms_colors('1234'))
      specrec = ef1D3_orig.specrec('10^3 km/s'); 
      flim = 10^-7;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end
    if 0 % fred mms4
      hca = irf_panel('reduced fe mms4');
      set(hca,'ColorOrder',mms_colors('1234'))
      specrec = ef1D4_orig.specrec('10^3 km/s'); 
      flim = 10^-7;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end
    if numel(ic_orig) == 3
      if 1 % B
        hca = irf_panel('Bx');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseB1.x,gseB2.x,gseB3.x},'comp');
        hca.YLabel.String = {'B_x','(nT)'};
        irf_legend(hca,{'1','2','3'},[0.02 0.98])
      end
      if 1 % n
        hca = irf_panel('n');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{ne1/n_lobe(1),ne2/n_lobe(2),ne3/n_lobe(3)},'comp')
        hca.YLabel.String = {'n/n^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])
        irf_legend(hca,{sprintf('n^{lb} = [%.3f,%.3f,%.3f] cc',n_lobe(1),n_lobe(2),n_lobe(3))},[0.02 0.98],'k')
      end
      if 1 % v par
        hca = irf_panel('v');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseVe1par/vte_lobe(1),gseVe2par/vte_lobe(2),gseVe3par/vte_lobe(3)},'comp');
        hca.YLabel.String = {'v_{e,||}/v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

        irf_legend(hca,{sprintf('v_{te}^{lb} = [%.0f,%.0f,%.0f] km/s',vte_lobe(1),vte_lobe(2),vte_lobe(3))},[0.02 0.98],'k')
      end
      if 1 % nv par
        hca = irf_panel('nv');
        set(hca,'ColorOrder',mms_colors('12345'))
        irf_plot(hca,{gseFlux1par/n_lobe(1)/vte_lobe(1),gseFlux2par/n_lobe(2)/vte_lobe(2),gseFlux3par/n_lobe(3)/vte_lobe(3)},'comp');
        hca.YLabel.String = {'nv_{||}/n^{lb}v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

      end
    else

      if 1 % B
        hca = irf_panel('Bx');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseB1.x,gseB2.x,gseB3.x,gseB4.x},'comp');
        hca.YLabel.String = {'B','(nT)'};
        irf_legend(hca,{'1','2','3','4'},[0.02 0.98])
      end
      if 1 % n
        hca = irf_panel('n');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{ne1/n_lobe(1),ne2/n_lobe(2),ne3/n_lobe(3),ne4/n_lobe(4)},'comp')
        hca.YLabel.String = {'n/n^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])
        irf_legend(hca,{sprintf('n^{lb} = [%.3f,%.3f,%.3f,%.3f] cc',n_lobe(1),n_lobe(2),n_lobe(3),n_lobe(4))},[0.02 0.98],'k')
      end
      if 1 % v par
        hca = irf_panel('v');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseVe1par/vte_lobe(1),gseVe2par/vte_lobe(2),gseVe3par/vte_lobe(3),gseVe4par/vte_lobe(4)},'comp');
        hca.YLabel.String = {'v_{e,||}/v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

        irf_legend(hca,{sprintf('v_{te}^{lb} = [%.0f,%.0f,%.0f,%.0f] km/s',vte_lobe(1),vte_lobe(2),vte_lobe(3),vte_lobe(4))},[0.02 0.98],'k')
      end
      if 1 % nv par
        hca = irf_panel('nv');
        set(hca,'ColorOrder',mms_colors('12345'))
        irf_plot(hca,{gseFlux1par/n_lobe(1)/vte_lobe(1),gseFlux2par/n_lobe(2)/vte_lobe(2),gseFlux3par/n_lobe(3)/vte_lobe(3),gseFlux4par/n_lobe(4)/vte_lobe(4)},'comp');
        hca.YLabel.String = {'nv_{||}/n^{lb}v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

      end
    end

    c_eval('irf_pl_mark(h(?),tint_phi,''y'')',2:nPanels)
    c_eval('irf_pl_mark(h(?),tint_sheet,''r'')',2:nPanels)
    c_eval('irf_pl_mark(h(?),tint_lobe,''b'')',2:nPanels)
    irf_zoom(h,'x',tint_fred)
    irf_zoom(h,'y')
    irf_plot_axis_align([h h_add])

    for ii = 1:numel(h_add)    
      h_add(ii).YLabel.Position(1) = 1.1;
    end

    cn.print(sprintf('tint%s-%s_flux_mms1234_map',tint_sep_utc(1,:),tint_sep_utc(2,:)),'path',printAccPotPath)
  end
  if 0 % all three together, mapping
    %%
    nPanels = 7;
    nRows = 3;
    nCols = 2;
    [h,h2] = initialize_combined_plot(nPanels,nRows,nCols,0.6,'vertical'); 
    
    linewidth = 2;
    if 1 % fred mms1
      hca = irf_panel('reduced fe mms1');
      set(hca,'ColorOrder',mms_colors('1234'))
      specrec = ef1D1_orig.specrec('10^3 km/s'); 
      flim = 10^-7;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end
    if numel(ic_orig) == 3
      if 1 % B
        hca = irf_panel('Bx');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseB1.x,gseB2.x,gseB3.x},'comp');
        hca.YLabel.String = {'B_x','(nT)'};
        irf_legend(hca,{'1','2','3'},[0.02 0.98])
      end
      if 1 % n
        hca = irf_panel('n');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{ne1/n_lobe(1),ne2/n_lobe(2),ne3/n_lobe(3)},'comp')
        hca.YLabel.String = {'n/n^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])
        irf_legend(hca,{sprintf('n^{lb} = [%.3f,%.3f,%.3f] cc',n_lobe(1),n_lobe(2),n_lobe(3))},[0.02 0.98],'k')
        hca.YTick = 0:1:20;
      end
      if 1 % v par
        hca = irf_panel('v');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseVe1par/vte_lobe(1),gseVe2par/vte_lobe(2),gseVe3par/vte_lobe(3)},'comp');
        hca.YLabel.String = {'v_{e,||}/v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

        irf_legend(hca,{sprintf('v_{te}^{lb} = [%.0f,%.0f,%.0f] km/s',vte_lobe(1),vte_lobe(2),vte_lobe(3))},[0.02 0.98],'k')
      end
      if 1 % nv par
        hca = irf_panel('nv');
        set(hca,'ColorOrder',mms_colors('12345'))
        irf_plot(hca,{gseFlux1par/n_lobe(1)/vte_lobe(1),gseFlux2par/n_lobe(2)/vte_lobe(2),gseFlux3par/n_lobe(3)/vte_lobe(3)},'comp');
        hca.YLabel.String = {'nv_{||}/n^{lb}v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

      end
      if 1 % map vtr
        hca = irf_panel('map vtr');
        set(hca,'ColorOrder',mms_colors('12345'))
        irf_plot(hca,{ts_vtr1,ts_vtr2,ts_vtr3},'comp');
        hold(hca,'on')
        %hlines = irf_plot(hca,{ts_vtr1flow_2,ts_vtr2flow_2,ts_vtr3flow_2},'comp');
        %c_eval('h_lines(?).LineWidth = linewidth;',1:numel(hlines))
        hold(hca,'off')
        hca.YLabel.String = {'(\psi/T_{e}^{lb})^{1/2}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])
      end
      if 1 % map vpsi
        hca = irf_panel('map vpsi');
        set(hca,'ColorOrder',mms_colors('12345'))
        irf_plot(hca,{ts_vpsi1,ts_vpsi2,ts_vpsi3},'comp');
        hold(hca,'on')
        %hlines = irf_plot(hca,{ts_vpsi1flow_2,ts_vpsi2flow_2,ts_vpsi3flow_2},'comp');
        %c_eval('h_lines(?).LineWidth = linewidth;',1:numel(hlines))
        hold(hca,'off')
        hca.YLabel.String = {'v_{\psi}/v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])
      end
    else
      if 1 % B
        hca = irf_panel('Bx');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseB1.x,gseB2.x,gseB3.x,gseB4.x},'comp');
        hca.YLabel.String = {'B','(nT)'};
        irf_legend(hca,{'1','2','3','4'},[0.02 0.98])
      end
      if 1 % n
        hca = irf_panel('n');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{ne1/n_lobe(1),ne2/n_lobe(2),ne3/n_lobe(3),ne4/n_lobe(4)},'comp')
        hca.YLabel.String = {'n/n^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])
        irf_legend(hca,{sprintf('n^{lb} = [%.3f,%.3f,%.3f,%.3f] cc',n_lobe(1),n_lobe(2),n_lobe(3),n_lobe(4))},[0.02 0.98],'k')
        hca.YTick = 0:1:20;
      end
      if 1 % v par
        hca = irf_panel('v');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseVe1par/vte_lobe(1),gseVe2par/vte_lobe(2),gseVe3par/vte_lobe(3),gseVe4par/vte_lobe(4)},'comp');
        hca.YLabel.String = {'v_{e,||}/v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

        irf_legend(hca,{sprintf('v_{te}^{lb} = [%.0f,%.0f,%.0f,%.0f] km/s',vte_lobe(1),vte_lobe(2),vte_lobe(3),vte_lobe(4))},[0.02 0.98],'k')
      end
      if 1 % nv par
        hca = irf_panel('nv');
        set(hca,'ColorOrder',mms_colors('12345'))
        irf_plot(hca,{gseFlux1par/n_lobe(1)/vte_lobe(1),gseFlux2par/n_lobe(2)/vte_lobe(2),gseFlux3par/n_lobe(3)/vte_lobe(3),gseFlux4par/n_lobe(4)/vte_lobe(4)},'comp');
        hca.YLabel.String = {'nv_{||}/n^{lb}v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

      end
      if 1 % map vtr
        hca = irf_panel('map vtr');
        set(hca,'ColorOrder',mms_colors('12345'))
        irf_plot(hca,{ts_vtr1,ts_vtr2,ts_vtr3,ts_vtr4},'comp');
        hold(hca,'on')
        %hlines = irf_plot(hca,{ts_vtr1flow_2,ts_vtr2flow_2,ts_vtr3flow_2,ts_vtr4flow_2},'comp');
        %c_eval('h_lines(?).LineWidth = linewidth;',1:numel(hlines))
        hold(hca,'off')        
        hca.YLabel.String = {'(\psi/T_{e}^{lb})^{1/2}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])
      end
      if 1 % map vpsi
        hca = irf_panel('map vpsi');
        set(hca,'ColorOrder',mms_colors('12345'))
        irf_plot(hca,{ts_vpsi1,ts_vpsi2,ts_vpsi3,ts_vpsi4},'comp');
        hold(hca,'on')
        %hlines = irf_plot(hca,{ts_vpsi1flow_2,ts_vpsi2flow_2,ts_vpsi3flow_2,ts_vpsi4flow_2},'comp');
        %c_eval('h_lines(?).LineWidth = linewidth;',1:numel(hlines))
        hold(hca,'off')        
        hca.YLabel.String = {'v_{\psi}/v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])
      end
    end

    for ii = 1:3
      hca = h2(ii);
      [h__,hcc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),Nmap,0.1:0.1:1,'k');
      clabel(h__,hcc);
      hold(hca,'on')
      [h__,hcc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),NVmap,0.2:0.2:3,'r');      
      clabel(h__,hcc);
      c_eval('plot(hca,ts_vtr?.data,ts_vpsi?.data,''ko'');',ii)
      c_eval('plot(hca,ts_vtr?flow.data,ts_vpsi?flow.data,''g*'');',ii)
      c_eval('plot(hca,ts_vtr?sep.data,ts_vpsi?sep.data,''y*'');',ii)
      hold(hca,'off')
      hca.XLabel.String = '(\psi/T_{e}^{lb})^{1/2}';
      hca.YLabel.String = 'v_{\psi}/v_{te}^{lb}';
      hca.Title.String = sprintf('mms %g',ii);
    end

    for ii = 1:3
      hca = h2(ii+3);
      [h__,hcc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),Nmap,0.1:0.1:1,'k');
      clabel(h__,hcc);
      hold(hca,'on')
      [h__,hcc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),NVmap,0.2:0.2:3,'r');      
      clabel(h__,hcc);
      %c_eval('plot(hca,ts_vtr?.data,ts_vpsi?.data,''ko'');',ii)
      c_eval('plot(hca,ts_vtr?flow_2.data,ts_vpsi?flow_2.data,''ko'');',ii)
      c_eval('plot(hca,ts_vtr?flow_2.data,ts_vpsi?flow_2.data,''g*'');',ii)
      c_eval('plot(hca,ts_vtr?sep_2.data,ts_vpsi?sep_2.data,''ko'');',ii)
      c_eval('plot(hca,ts_vtr?sep_2.data,ts_vpsi?sep_2.data,''y*'');',ii)
      hold(hca,'off')
      hca.XLabel.String = '(\psi/T_{e}^{lb})^{1/2}';
      hca.YLabel.String = 'v_{\psi}/v_{te}^{lb}';
      hca.Title.String = sprintf('mms %g',ii);
    end
    h2(4).Title.String = {'flux mult by 2',h2(4).Title.String};
    
    c_eval('irf_pl_mark(h(?),tint_sep,''y'')',2:nPanels)
    c_eval('irf_pl_mark(h(?),tint_flow,''g'')',2:nPanels)
    c_eval('irf_pl_mark(h(?),tint_sheet,''r'')',2:nPanels)
    c_eval('irf_pl_mark(h(?),tint_lobe,''b'')',2:nPanels)
    irf_zoom(h,'x',tint_fred)
    irf_zoom(h,'y')
    irf_plot_axis_align(h)
    

    cn.print(sprintf('tint%s-%s_flux_mms1234_map_event%g',tint_sep_utc(1,:),tint_sep_utc(2,:),events(ievent)),'path',[printAccPotPath 'flux1234_map/'])    
  end
  if 0 % Make figures, all three/four together  
    nPanels = 5;
    h = irf_plot(nPanels); 
    h_add = h;
    h_add(:) = [];
    if 1 % fred mms1
      hca = irf_panel('reduced fe mms1');
      set(hca,'ColorOrder',mms_colors('1234'))
      specrec = ef1D1_orig.specrec('10^3 km/s'); 
      flim = 10^-7;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end
    if 0 % fred mms2
      hca = irf_panel('reduced fe mms2');
      set(hca,'ColorOrder',mms_colors('1234'))
      specrec = ef1D2_orig.specrec('10^3 km/s'); 
      flim = 10^-7;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end
    if 0 % fred mms1
      hca = irf_panel('reduced fe mms3');
      set(hca,'ColorOrder',mms_colors('1234'))
      specrec = ef1D3_orig.specrec('10^3 km/s'); 
      flim = 10^-7;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end
    if 0 % fred mms4
      hca = irf_panel('reduced fe mms4');
      set(hca,'ColorOrder',mms_colors('1234'))
      specrec = ef1D4_orig.specrec('10^3 km/s'); 
      flim = 10^-7;
      specrec.p(specrec.p<flim) = NaN;
      irf_spectrogram(hca,specrec);
      hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};
      %irf_legend(hca,legends_mms,[0.02 0.98])
    end
    if numel(ic_orig) == 3
      if 1 % B
        hca = irf_panel('Bx');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseB1.x,gseB2.x,gseB3.x},'comp');
        hca.YLabel.String = {'B_x','(nT)'};
        irf_legend(hca,{'1','2','3'},[0.02 0.98])
      end
      if 1 % n
        hca = irf_panel('n');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{ne1/n_lobe(1),ne2/n_lobe(2),ne3/n_lobe(3)},'comp')
        hca.YLabel.String = {'n/n^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])
        irf_legend(hca,{sprintf('n^{lb} = [%.3f,%.3f,%.3f] cc',n_lobe(1),n_lobe(2),n_lobe(3))},[0.02 0.98],'k')
      end
      if 1 % v par
        hca = irf_panel('v');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseVe1par/vte_lobe(1),gseVe2par/vte_lobe(2),gseVe3par/vte_lobe(3)},'comp');
        hca.YLabel.String = {'v_{e,||}/v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

        irf_legend(hca,{sprintf('v_{te}^{lb} = [%.0f,%.0f,%.0f] km/s',vte_lobe(1),vte_lobe(2),vte_lobe(3))},[0.02 0.98],'k')
      end
      if 1 % nv par
        hca = irf_panel('nv');
        set(hca,'ColorOrder',mms_colors('12345'))
        irf_plot(hca,{gseFlux1par/n_lobe(1)/vte_lobe(1),gseFlux2par/n_lobe(2)/vte_lobe(2),gseFlux3par/n_lobe(3)/vte_lobe(3)},'comp');
        hca.YLabel.String = {'nv_{||}/n^{lb}v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

      end
    else

      if 1 % B
        hca = irf_panel('Bx');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseB1.x,gseB2.x,gseB3.x,gseB4.x},'comp');
        hca.YLabel.String = {'B','(nT)'};
        irf_legend(hca,{'1','2','3','4'},[0.02 0.98])
      end
      if 1 % n
        hca = irf_panel('n');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{ne1/n_lobe(1),ne2/n_lobe(2),ne3/n_lobe(3),ne4/n_lobe(4)},'comp')
        hca.YLabel.String = {'n/n^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])
        irf_legend(hca,{sprintf('n^{lb} = [%.3f,%.3f,%.3f,%.3f] cc',n_lobe(1),n_lobe(2),n_lobe(3),n_lobe(4))},[0.02 0.98],'k')
      end
      if 1 % v par
        hca = irf_panel('v');
        set(hca,'ColorOrder',mms_colors('1234'))
        irf_plot(hca,{gseVe1par/vte_lobe(1),gseVe2par/vte_lobe(2),gseVe3par/vte_lobe(3),gseVe4par/vte_lobe(4)},'comp');
        hca.YLabel.String = {'v_{e,||}/v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

        irf_legend(hca,{sprintf('v_{te}^{lb} = [%.0f,%.0f,%.0f,%.0f] km/s',vte_lobe(1),vte_lobe(2),vte_lobe(3),vte_lobe(4))},[0.02 0.98],'k')
      end
      if 1 % nv par
        hca = irf_panel('nv');
        set(hca,'ColorOrder',mms_colors('12345'))
        irf_plot(hca,{gseFlux1par/n_lobe(1)/vte_lobe(1),gseFlux2par/n_lobe(2)/vte_lobe(2),gseFlux3par/n_lobe(3)/vte_lobe(3),gseFlux4par/n_lobe(4)/vte_lobe(4)},'comp');
        hca.YLabel.String = {'nv_{||}/n^{lb}v_{te}^{lb}',''};
        hca.YLabel.Interpreter = 'tex';
        %irf_legend(hca,{'x','y','z'},[0.02 0.98])

      end
    end

    c_eval('irf_pl_mark(h(?),tint_phi,''y'')',2:nPanels)
    c_eval('irf_pl_mark(h(?),tint_sheet,''r'')',2:nPanels)
    c_eval('irf_pl_mark(h(?),tint_lobe,''b'')',2:nPanels)
    irf_zoom(h,'x',tint_fred)
    irf_zoom(h,'y')
    irf_plot_axis_align([h h_add])

    for ii = 1:numel(h_add)    
      h_add(ii).YLabel.Position(1) = 1.1;
    end

    cn.print(sprintf('tint%s-%s_flux_mms1234',tint_sep_utc(1,:),tint_sep_utc(2,:)),'path',printAccPotPath)
  end
  %% Save data to file
  if 0
  % .mat or .txt (together or separate?)  
  acc_nv_data.event_id = event;
  acc_nv_data.ic = ic;  
  acc_nv_data.n_lobe = n_lobe;
  acc_nv_data.n_sheet = n_sheet;
  acc_nv_data.n_sep = n_sep;
  acc_nv_data.n_sep_min = n_sep_min;
  
  save(sprintf('%s/acc_nv_data_event_%g',saveAccPotPath,event),'acc_nv_data')
  %fid = fopen([matlabPath 'acc_pot_statistics.txt'],'a+');
  % event_id
  %save_format = '%f';
  %fid = fclose(fid);
  %fprintf('printed to file %sesw_properties.txt: %s /n',matlabPath,print_str)
  end
  %close all
  clear tint_flow
end