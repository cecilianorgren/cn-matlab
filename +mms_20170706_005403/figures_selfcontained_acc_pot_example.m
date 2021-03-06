units = irf_units;

%% Time intervals
% burst interval
tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');
% shorter time interval for dispersion analysis
Tints = irf.tint('2017-07-06T00:55:39.50Z/2017-07-06T00:55:41.50Z'); % second "good" batch
Tints = irf.tint('2017-07-06T00:54:13.50Z/2017-07-06T00:54:17.00Z'); % first "good" batch
Tints = irf.tint('2017-07-06T00:54:13.50Z/2017-07-06T00:54:17.00Z'); % acceleration channel, first "good" batch
event = 3;
sep.get_tints;

%% Load data
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)
c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
%c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)
R = mms.get_data('R_gse',tint);
if size(R.gseR1,2) == 4
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
else
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
end

c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)

c_eval('tic; ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?); toc;',ic);
c_eval('tic; ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?); toc;',ic);

c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)

c_eval('beta?e = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
c_eval('beta?i = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
c_eval('beta? = beta?i + beta?e.resample(beta?i);',ic)

%% Pick out parameters in the lobe, sheet, and separatrix intervals
c_eval('B_lobe = mean(gseB?.abs.tlim(tint_lobe).data,1);',ic)

c_eval('n_lobe = mean(ne?.tlim(tint_lobe).data,1);',ic)
c_eval('n_sheet = mean(ne?.tlim(tint_sheet).data,1);',ic)
c_eval('n_sep = mean(ne?.tlim(tint_sep).data,1);',ic)

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

%% % Load wave phase velocities (obtained manually and semi-manually)
fid = fopen([matlabPath 'esw_properties_redo.txt'],'r');
data_format_read = '%s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f'; 
esw_data = textscan(fid,[data_format_read]);
fclose(fid)

if 1 % load esw_properties_irf_4_v_gui.txt
  fid = fopen([matlabPath 'esw_properties_irf_4_v_gui.txt'],'r');
  data_format_read = '%s %s %s %s %s %f %f %f %f %f %f %f %f %f %f %f %f %f'; 
  esw_data_ = textscan(fid,[data_format_read]);
  fclose(fid)
  for icell = 1:numel(esw_data_)
    esw_data{icell} = cat(1,esw_data{icell},esw_data_{icell});
  end
end

all_t_center = EpochTT(char(esw_data{5}));
[all_t_center_sorted,ind_sorted] = all_t_center.sort;

for icell = 1:numel(esw_data)
  esw_data{icell} = esw_data{icell}(ind_sorted);  
end

tsVph = irf.ts_vec_xyz(char(esw_data{5}),[esw_data{6} esw_data{7} esw_data{8}]);
tsVphpar = tsVph.dot(gseB1.norm.resample(tsVph));
tsPhi = irf.ts_scalar(char(esw_data{5}),[esw_data{10} esw_data{11} esw_data{12} esw_data{13}]);
c_eval('tsVtrap? = irf.ts_scalar(char(esw_data{5}),sqrt(2*units.e*esw_data{9+?}/units.me)*1e-3);',1:4)
c_eval('vmin? = tsVphpar-tsVtrap?;',1:4)
c_eval('vmax? = tsVphpar+tsVtrap?;',1:4)
vmax = 0.25*(vmax1 + vmax2 + vmax3 + vmax4); vmax.name = 'av(vph+vtrap)';
vmin = 0.25*(vmin1 + vmin2 + vmin3 + vmin4); vmin.name = 'av(vph-vtrap)';
  
%% Reduced electron distribution
eint = [00 40000];
lowerelim = 000;
vint = [-Inf Inf];
tint_fred = tint_fred;%tint_phi;
c_eval('eDist = ePDist?.tlim(tint_fred);',ic)
eDist_orig = eDist;

% remove background
nSecondary = 5;
nPhoto = 1;
eDist_nobg = mms.remove_edist_background(eDist,'nSecondary',nSecondary,'Nphotoe_art',nPhoto);
%eDist = eDist_nobg;

c_eval('scpot = scPot?.resample(eDist);',ic)
c_eval('dir_red = dmpaB?.resample(eDist).norm;',ic)
energies = sqrt(2*eDist.depend{1}(1,:)*units.eV/units.me)/1000; % km/s
vgmax = 70000;
vg = -vgmax:1000:vgmax;
vg(abs(vg)>70000) = [];

nMC = 500;
%tic; ef1D = eDist.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'vg',vg); toc % reduced distribution along B
tic; ef1D_orig = eDist_orig.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'nMC',nMC,'vg',vg); toc % reduced distribution along B
tic; ef1D_nobg = eDist_nobg.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'nMC',nMC,'vg',vg); toc % reduced distribution along B
%tic; ef1D_new = eDist_new.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'nMC',nMC,'vg',vg); toc % reduced distribution along B

%ef1D_orig = ef1D;
%ef1D = ef1D_new;
if 0
  %%
  fsheet_orig = ef1D_orig.tlim(tint_sheet);
  fsheet_nobg = ef1D_nobg.tlim(tint_sheet);
  plot(fsheet_orig.depend{1}(1,:),mean(fsheet_orig.data,1),...
       fsheet_nobg.depend{1}(1,:),mean(fsheet_nobg.data,1))
end

%% Plot acceleration potential as a function of energy
fphi_orig = ef1D_orig.tlim(tint_phi);
fphi_nobg = ef1D_nobg.tlim(tint_phi);

fphi = fphi_nobg;

energies_red = sign(ef1D.depend{1}(1,:)).*(ef1D.depend{1}(1,:)*1000).^2*units.me/2/units.eV; % km/s
plot_x = sign(energies_red).*log10(energies_red);
plot_x = energies_red*1e-3;

hca = subplot(2,1,1);
plot(hca,plot_x,fphi.data(25:1:end-10,:)')
%hca.XLim = 0.05*[-1 1];
%hca.YLim = [0 1]*1e-3;
%hca.XScale = 'log';
hca = subplot(2,1,2);
plot(hca,fphi.depend{1}(1,:)*1e-3,fphi.data(10:1:end-10,:)')
hca.XLabel.String = 'v (10^3 km/s)';
%hca.XLim = 10*[-1 1];
%hca.YLim = [0 1]*1e-3;

%% Plot acceleration potential and electron particle distributions
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?.resample(gseVe?);',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
%c_eval('ePitch = ePitch?;',ic)

fred = ef1D_nobg;

h2 = setup_subplots(6,1);

times_phi_utc = fred.tlim(tint_phi).time.utc; times_phi_utc = times_phi_utc(:,12:23);

fontsize = 8;
  
colors = mms_colors('xzy');
isub = 1;
if 1 % average fred during sep, lobe and sheet interval
  hca = h2(isub); isub = isub + 1;
  
  if 1
    fered_sep = fred.tlim(tint_sep);
    plot_v_sep = nanmean(fered_sep.depend{1},1);
    plot_f_sep = nanmean(fered_sep.data,1);
  end
  fered_lobe = fred.tlim(tint_lobe);
  plot_v_lobe = nanmean(fered_lobe.depend{1},1);
  plot_f_lobe = nanmean(fered_lobe.data,1);
  
  fered_sheet = fred.tlim(tint_sheet);
  plot_v_sheet = nanmean(fered_sheet.depend{1},1);
  plot_f_sheet = nanmean(fered_sheet.data,1);
  
  h_fered = plot(hca,plot_v_lobe*1e-3,plot_f_lobe,...
                     plot_v_sheet*1e-3,plot_f_sheet);
                   
  c_eval('h_fered(?).Color = colors(?,:);',1:2)
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = sprintf('f_e (%s)',fered_sep.units);
  hca.XLim = [-50 50];
  hca.XTick = [-60:20:60];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  set(hca,'colororder',colors);
  tint_str = {sprintf('lobe: %s - %s',tint_lobe_utc(1,:),tint_lobe_utc(2,:));...
              sprintf('plasma sheet: %s - %s',tint_sheet_utc(1,:),tint_sheet_utc(2,:))};
  irf_legend(hca,tint_str,[0.01 0.98],'fontsize',fontsize)
  hca.Title.String = sprintf('E > %g eV',lowerelim);
end

%ndist_to_combine = 3;
n_int = 5;
lowerelim_peaks = 400;
Elim_axis = 4000;
colors = mms_colors('matlab');
colors = [colors(3:end,:); colors(1:2,:)];
clear plot_f_phi f_sep_utc hmark_all str_Emax Emax Emax_within_interval 
clear Emax_av Emax_av_ind plot_f_phi_Emax Emax_utc f_phi_utc

if 1 % divide f into intervals, common for both v,E x-axis, mark TSeries
  fred_phi = fred.tlim(tint_phi);
  plot_v_phi = mean(fred_phi.depend{1},1);
  plot_E_phi = sign(plot_v_phi)*0.5*units.me.*(plot_v_phi*1e3).^2/units.e;   
  
  fred_phi_peaks = fred_phi;   
  fred_phi_peaks.data(abs(repmat(plot_E_phi,fred_phi_peaks.length,1))<lowerelim_peaks) = 0;
  [val,ind] = max(fred_phi_peaks.data');
  Emax = (plot_E_phi(ind));
     
  nf = fred_phi.length;
  [ind1,ind2] = ind12(nf,n_int);
  for i_int = 1:n_int
    f_phi_utc{i_int,1} = sprintf('%s - %s',times_phi_utc(ind1(i_int),:),times_phi_utc(ind2(i_int),:));
    Emax_tmp = Emax(ind1(i_int):ind2(i_int));
    [val_tmp_Emax,ind_tmp_Emax] = max(abs(Emax_tmp));
    Emax_within_interval(i_int) = Emax_tmp(ind_tmp_Emax);
    time_Emax_within_interval{i_int} = fred_phi_peaks.time(ind1(i_int)-1+ind_tmp_Emax);
    %Emax_utc = 
        
    plot_f_phi(i_int,:) = mean(fred_phi.data(ind1(i_int):ind2(i_int),:),1);    
    plot_f_phi_Emax(i_int,:) = fred_phi.data(ind1(i_int)+ind_tmp_Emax-1,:);
    
    Emax_tmp = plot_f_phi(i_int,:); Emax_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(Emax_tmp);
    Emax_av(i_int) = (plot_E_phi(ind));
    Emax_av_ind(i_int) = ind;
    f_at_Emax_av(i_int) = plot_f_phi(i_int,ind);
    
    str_Emax{i_int,1} = sprintf('phimax =[%.0f, %.0f] eV',Emax_av(i_int),Emax_within_interval(i_int));
    
    if 0
    [hmark,cc] = irf_pl_mark(h(npanels_large + 2),fred_phi.time([ind1(i_int) ind2(i_int)]).epochUnix');    
    if isa(hmark,'matlab.graphics.primitive.Line')
      hmark.Color = colors(i_int,:);
    elseif isa(hmark,'matlab.graphics.primitive.Patch')
      hmark.YData = hmark.YData + 120;
      hmark.ZData = [1 1 1 1];
      hmark.EdgeColor = [0 0 0];
      hmark.FaceColor = colors(i_int,:);
    end    
    hmark_all{i_int} = hmark;
    end
  end  
end

if 1 % all fred during sep interval, vs vpar
  hca = h2(isub); isub = isub + 1;
  h_fered = plot(hca,plot_v_sep*1e-3,plot_f_phi);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)  
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = sprintf('f_e (%s)',fered_sep.units);
  hca.XLim = [-30 50];
  set(hca,'colororder',colors)  
  hca.Title.String = 'Average over indicated tint';
  if plot_E_phi(ind) > 0
    irf_legend(hca,f_phi_utc,[0.01 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,f_phi_utc,[0.98 0.98],'fontsize',fontsize)
  end
end
if 1 % all fred during sep interval, vs +-Epar
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_E_phi*1e0,plot_f_phi);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    h_fered_star = plot(hca,Emax_av(i_int),f_at_Emax_av(i_int),'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'E_{||} (eV)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  hca.XLim = Elim_axis*[-1 1]+ 0.5*mean(Emax_av);
  hca.XTick = -30000:1000:30000;
     
  hca.Title.String = {'Average over indicated tint',sprintf('phimax with lower limit E > %g eV',lowerelim_peaks)};

  if 1 % plot lowerelim_peaks
    hold(hca,'on')
    plot(hca,lowerelim_peaks*[-1 -1],hca.YLim,'--k',lowerelim_peaks*[1 1],hca.YLim,'--k')
    hold(hca,'off')
  end
  
  set(hca,'colororder',colors)
  if plot_E_phi(ind) > 0
    irf_legend(hca,str_Emax,[0.02 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,str_Emax,[0.98 0.98],'fontsize',fontsize)
  end
end
if 1 % all fred during sep interval, vs vpar, max within interal
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_v_sep,plot_f_phi_Emax);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    f_tmp = plot_f_phi_Emax(i_int,:);
    f_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(abs(f_tmp));
    h_fered_star = plot(hca,plot_v_sep(ind),val,'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  %hca.XLim = Elim_axis*[-1 1]+ 0.5*mean(Emax);
  %hca.XTick = -30000:1000:30000;
  if hca.XLim(2)<max(Emax)
    hca.XLim = max(Emax)*1*[-1 1] + 0.5*max(Emax);
    hca.XTick = -30000:2000:30000;
  end
     
  hca.Title.String = sprintf('E > %g eV',lowerelim_peaks);
  if 1 % plot lowerelim_peaks
    hold(hca,'on')
    plot(hca,lowerelim_peaks*[-1 -1],hca.YLim,'--k',lowerelim_peaks*[1 1],hca.YLim,'--k')
    hold(hca,'off')
  end
  
  set(hca,'colororder',colors)
  if plot_E_phi(ind) > 0
    irf_legend(hca,str_Emax,[0.02 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,str_Emax,[0.98 0.98],'fontsize',fontsize)
  end
  hca.Title.String = {'Max phi in indicated tint',sprintf('phimax with lower limit E > %g eV',lowerelim_peaks)};
end
if 1 % all fred during sep interval, vs +-Epar, max within interal
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_E_phi,plot_f_phi_Emax);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    f_tmp = plot_f_phi_Emax(i_int,:);
    f_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(abs(f_tmp));
    h_fered_star = plot(hca,plot_E_phi(ind),val,'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'E_{||} (eV)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  hca.XLim = Elim_axis*[-1 1]+ 0.5*mean(Emax);
  hca.XTick = -30000:1000:30000;
  if hca.XLim(2)<max(Emax)
    hca.XLim = max(Emax)*1*[-1 1] + 0.5*max(Emax);
    hca.XTick = -30000:2000:30000;
  end
  hca.XLim = [-2000 8000];
     
  hca.Title.String = sprintf('E > %g eV',lowerelim_peaks);
  if 1 % plot lowerelim_peaks
    hold(hca,'on')
    plot(hca,lowerelim_peaks*[-1 -1],hca.YLim,'--k',lowerelim_peaks*[1 1],hca.YLim,'--k')
    hold(hca,'off')
  end
  
  set(hca,'colororder',colors)
  if plot_E_phi(ind) > 0
    irf_legend(hca,str_Emax,[0.02 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,str_Emax,[0.98 0.98],'fontsize',fontsize)
  end
  hca.Title.String = {'Max phi in indicated tint',sprintf('phimax with lower limit E > %g eV',lowerelim_peaks)};
end

if 1 % all fred during sep interval, vs +-Epar, max within interal
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_E_phi,plot_f_phi_Emax);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    f_tmp = plot_f_phi_Emax(i_int,:);
    f_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(abs(f_tmp));
    h_fered_star = plot(hca,plot_E_phi(ind),val,'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'E_{||} (eV)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  hca.XLim = Elim_axis*[-1 1]+ 0.5*mean(Emax);
  hca.XTick = -30000:1000:30000;
  if hca.XLim(2)<max(Emax)
    hca.XLim = max(Emax)*1*[-1 1] + 0.5*max(Emax);
    hca.XTick = -30000:2000:30000;
  end
     
  hca.Title.String = sprintf('E > %g eV',lowerelim_peaks);
  if 1 % plot lowerelim_peaks
    hold(hca,'on')
    plot(hca,lowerelim_peaks*[-1 -1],hca.YLim,'--k',lowerelim_peaks*[1 1],hca.YLim,'--k')
    hold(hca,'off')
  end
  
  set(hca,'colororder',colors)
  if 0
    if plot_E_phi(ind) > 0
      irf_legend(hca,str_Emax,[0.02 0.98],'fontsize',fontsize)
    else
      irf_legend(hca,str_Emax,[0.98 0.98],'fontsize',fontsize)
    end
  end
  legs_Emax = {};
  for i_int = 1:n_int
    time_utc_tmp = time_Emax_within_interval{i_int}.utc;
    %legs_Emax{i_int,1} = sprintf('%s',time_utc_tmp(12:23));
    legs_Emax{i_int,1} = sprintf('%s: acc. pot. = %g eV',time_utc_tmp(12:23),round(Emax_within_interval(i_int)/100)*100);
  end
  irf_legend(hca,legs_Emax,[0.98 0.98],'fontsize',fontsize)
  %hca.Title.String = {'Max phi in indicated tint',sprintf('phimax with lower limit E > %g eV',lowerelim_peaks)};
  hca.Title.String = '';
  hca.YLim = [0 0.6]*1e-3;
  hca.XLim = [-2000 6000];  
end

if 0 % all lobe, plasma sheet, fred during sep interval, vs vpar
  hca = h2(isub); isub = isub + 1;
  if 1 % acceleration channel, 
    h_fered = plot(hca,plot_v_sep*1e-3,plot_f_phi);
    c_eval('h_fered(?).Color = colors(?,:);',1:n_int)  
  end
  if 1
    hold(hca,'on')
    fered_lobe = ef1D.tlim(tint_lobe);
    plot_v_lobe = nanmean(fered_lobe.depend{1},1);
    plot_f_lobe = nanmean(fered_lobe.data,1);

    fered_sheet = ef1D.tlim(tint_sheet);
    plot_v_sheet = nanmean(fered_sheet.depend{1},1);
    plot_f_sheet = nanmean(fered_sheet.data,1);

    h_fered = plot(hca,plot_v_lobe*1e-3,plot_f_lobe/5,...
                       plot_v_sheet*1e-3,plot_f_sheet/15);

    
    h_fered(1).Color = mms_colors('x');
    h_fered(2).Color = mms_colors('z');
    hold(hca,'off')
  end
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = sprintf('f_e (%s)',fered_sep.units);
  hca.XLim = [-40 40];
  set(hca,'colororder',colors)  
  hca.Title.String = 'Average over indicated tint';
  if plot_E_phi(ind) > 0
    irf_legend(hca,f_phi_utc,[0.01 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,f_phi_utc,[0.98 0.98],'fontsize',fontsize)
  end
end

n_intervals_sep = 5;
vgmax = 60000;
vg = -vgmax:2000:vgmax;
%lowerelim = 40;
if 0 % 1 reduced distribution, 1st half of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[0 1];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  
  c_eval('ePara = dmpaB?slow.resample(eDist).norm;',ic)
  c_eval('ePerp1 = ePara.cross(dslE?slow.resample(eDist)).norm;',ic)
  ePerp2 = ePara.cross(ePerp1).norm;
  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  %axis(hca,'equal')
  axis(hca,'square')
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
end
if 0 % 1 reduced distribution, 2nd half of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[1 2];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  

  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  %axis(hca,'equal')
  axis(hca,'square')
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
end
if 0 % 1 reduced distribution, 3rd third of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[2 3];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  
  c_eval('ePara = dmpaB?slow.resample(eDist).norm;',ic)
  c_eval('ePerp1 = ePara.cross(dslE?slow.resample(eDist)).norm;',ic)
  ePerp2 = ePara.cross(ePerp1).norm;
  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
  %axis(hca,'equal')
  axis(hca,'square')
end

for ih = 1:numel(h2)
  h2(ih).FontSize = fontsize;
end
%for ih = 4:6, h2(ih).CLim = [-15 -10.2]; end

%% Plot acceleration potential and electron particle distributions, only what is needed for paper
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?.resample(gseVe?);',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
%c_eval('ePitch = ePitch?;',ic)

n_int = 6;
lowerelim_peaks = 400;
Elim_axis = 4000;
colors = mms_colors('matlab');
colors = [colors(3:end,:); colors(1:2,:)];
clear plot_f_phi f_sep_utc hmark_all str_Emax Emax Emax_within_interval 
clear Emax_av Emax_av_ind plot_f_phi_Emax Emax_utc f_phi_utc

if 1 % divide f into intervals, common for both v,E x-axis, mark TSeries
  fred_phi = ef1D.tlim(tint_phi);
  plot_v_phi = mean(fred_phi.depend{1},1);
  plot_E_phi = sign(plot_v_phi)*0.5*units.me.*(plot_v_phi*1e3).^2/units.e;   
  
  fred_phi_peaks = fred_phi;   
  fred_phi_peaks.data(abs(repmat(plot_E_phi,fred_phi_peaks.length,1))<lowerelim_peaks) = 0;
  [val,ind] = max(fred_phi_peaks.data');
  Emax = (plot_E_phi(ind));
     
  nf = fred_phi.length;
  [ind1,ind2] = ind12(nf,n_int);
  for i_int = 1:n_int
    f_phi_utc{i_int,1} = sprintf('%s - %s',times_phi_utc(ind1(i_int),:),times_phi_utc(ind2(i_int),:));
    Emax_tmp = Emax(ind1(i_int):ind2(i_int));
    [val_tmp_Emax,ind_tmp_Emax] = max(abs(Emax_tmp));
    Emax_within_interval(i_int) = Emax_tmp(ind_tmp_Emax);
    time_Emax_within_interval{i_int} = fred_phi_peaks.time(ind1(i_int)-1+ind_tmp_Emax);
    %Emax_utc = 
        
    plot_f_phi(i_int,:) = mean(fred_phi.data(ind1(i_int):ind2(i_int),:),1);    
    plot_f_phi_Emax(i_int,:) = fred_phi.data(ind1(i_int)+ind_tmp_Emax-1,:);
    
    Emax_tmp = plot_f_phi(i_int,:); Emax_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(Emax_tmp);
    Emax_av(i_int) = (plot_E_phi(ind));
    Emax_av_ind(i_int) = ind;
    f_at_Emax_av(i_int) = plot_f_phi(i_int,ind);
    
    str_Emax{i_int,1} = sprintf('phimax =[%.0f, %.0f] eV',Emax_av(i_int),Emax_within_interval(i_int));
    
    if 0
    [hmark,cc] = irf_pl_mark(h(npanels_large + 2),fred_phi.time([ind1(i_int) ind2(i_int)]).epochUnix');    
    if isa(hmark,'matlab.graphics.primitive.Line')
      hmark.Color = colors(i_int,:);
    elseif isa(hmark,'matlab.graphics.primitive.Patch')
      hmark.YData = hmark.YData + 120;
      hmark.ZData = [1 1 1 1];
      hmark.EdgeColor = [0 0 0];
      hmark.FaceColor = colors(i_int,:);
    end    
    hmark_all{i_int} = hmark;
    end
  end  
end

h2 = setup_subplots(3,1);
%h2().Position(2) = 0.2;
%h2.Position(4) = 0.7;

times_phi_utc = ef1D.tlim(tint_phi).time.utc; times_phi_utc = times_phi_utc(:,12:23);

fontsize = 12;
  
colors = mms_colors('xzy');
isub = 1;
if 1 % time series to mark time interval and chosen phi's
  hca = h2(isub); isub = isub +1;  
  fred_min = 1e-6;
  fred_to_plot = ef1D.tlim(tint_phi); fred_to_plot.data(fred_to_plot.data < fred_min) = NaN;
  irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));  
  hca.YLim = fred_to_plot.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'}; 
  irf_timeaxis(hca);
end
if 0 % average fred during sep, lobe and sheet interval
  hca = h2(isub); isub = isub + 1;
  
  if 0
    fered_sep = ef1D.tlim(tint_sep);
    plot_v_sep = nanmean(fered_sep.depend{1},1);
    plot_f_sep = nanmean(fered_sep.data,1);
  end
  fered_lobe = ef1D.tlim(tint_lobe);
  plot_v_lobe = nanmean(fered_lobe.depend{1},1);
  plot_f_lobe = nanmean(fered_lobe.data,1);
  
  fered_sheet = ef1D.tlim(tint_sheet);
  plot_v_sheet = nanmean(fered_sheet.depend{1},1);
  plot_f_sheet = nanmean(fered_sheet.data,1);
  
  h_fered = plot(hca,plot_v_lobe*1e-3,plot_f_lobe,...
                     plot_v_sheet*1e-3,plot_f_sheet);
                   
  c_eval('h_fered(?).Color = colors(?,:);',1:2)
  hca.XLabel.String = 'v (10^3 km/s)';
  hca.YLabel.String = sprintf('f_e (%s)',fered_sep.units);
  hca.XLim = [-50 50];
  hca.XTick = [-60:20:60];
  
  set(hca,'colororder',colors);
  tint_str = {sprintf('lobe: %s - %s',tint_lobe_utc(1,:),tint_lobe_utc(2,:));...
              sprintf('plasma sheet: %s - %s',tint_sheet_utc(1,:),tint_sheet_utc(2,:))};
  irf_legend(hca,tint_str,[0.01 0.98],'fontsize',fontsize)
  hca.Title.String = sprintf('E > %g eV',lowerelim);
end

%ndist_to_combine = 3;


if 1 % all fred during sep interval, vs +-Epar, max within interal
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_E_phi,plot_f_phi);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    f_tmp = plot_f_phi_Emax(i_int,:);
    f_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(abs(f_tmp));
    h_fered_star = plot(hca,Emax_av(i_int),f_at_Emax_av(i_int),'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'E_{||} (eV)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  hca.XLim = Elim_axis*[-1 1]+ 0.5*mean(Emax);
  hca.XTick = -30000:1000:30000;
  if hca.XLim(2)<max(Emax)
    hca.XLim = max(Emax)*1*[-1 1] + 0.5*max(Emax);
    hca.XTick = -30000:2000:30000;
  end
     
  hca.Title.String = sprintf('E > %g eV',lowerelim_peaks);
  if 0 % plot lowerelim_peaks
    hold(hca,'on')
    plot(hca,lowerelim_peaks*[-1 -1],hca.YLim,'--k',lowerelim_peaks*[1 1],hca.YLim,'--k')
    hold(hca,'off')
  end
  
  set(hca,'colororder',colors)
  if 0
    if plot_E_phi(ind) > 0
      irf_legend(hca,str_Emax,[0.02 0.98],'fontsize',fontsize)
    else
      irf_legend(hca,str_Emax,[0.98 0.98],'fontsize',fontsize)
    end
  end
  legs_Emax = {};
  for i_int = 1:n_int
    time_utc_tmp = time_Emax_within_interval{i_int}.utc;
    %legs_Emax{i_int,1} = sprintf('%s',time_utc_tmp(12:23));
    %irf_legend(hca,,[0.01 0.98],'fontsize',fontsize)
    legs_Emax{i_int,1} = sprintf('%s: acc. pot. = %4.0f eV',f_phi_utc{i_int},round(Emax_av(i_int)/100)*100);
  end
  
  irf_legend(hca,legs_Emax,[0.98 0.98],'fontsize',fontsize)
  %hca.Title.String = {'Max phi in indicated tint',sprintf('phimax with lower limit E > %g eV',lowerelim_peaks)};
  hca.Title.String = '';
  hca.YLim = [0 0.6]*1e-3;
  hca.XLim = [-2000 6000];  
end
if 1 % all fred during sep interval, vs +-Epar, max within interal
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_E_phi,plot_f_phi_Emax);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    f_tmp = plot_f_phi_Emax(i_int,:);
    f_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(abs(f_tmp));
    h_fered_star = plot(hca,plot_E_phi(ind),val,'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'E_{||} (eV)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  hca.Title.String = '';           
  
  set(hca,'colororder',colors)  
  legs_Emax = {};
  for i_int = 1:n_int
    time_utc_tmp = time_Emax_within_interval{i_int}.utc;
    %legs_Emax{i_int,1} = sprintf('%s',time_utc_tmp(12:23));
    legs_Emax{i_int,1} = sprintf('%s: acc. pot. = %g eV',time_utc_tmp(12:23),round(Emax_within_interval(i_int)/100)*100);
  end
  irf_legend(hca,legs_Emax,[0.98 0.98],'fontsize',fontsize) 
  
  hca.YLim = [0 0.6]*1e-3;
  hca.XLim = [-2000 6000];
  hca.XTick = -30000:1000:30000;
end


n_intervals_sep = 5;
vgmax = 60000;
vg = -vgmax:2000:vgmax;
%lowerelim = 40;
if 0 % 1 reduced distribution, 1st half of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[0 1];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  
  c_eval('ePara = dmpaB?slow.resample(eDist).norm;',ic)
  c_eval('ePerp1 = ePara.cross(dslE?slow.resample(eDist)).norm;',ic)
  ePerp2 = ePara.cross(ePerp1).norm;
  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  %axis(hca,'equal')
  axis(hca,'square')
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
end
if 0 % 1 reduced distribution, 2nd half of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[1 2];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  

  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  %axis(hca,'equal')
  axis(hca,'square')
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
end
if 0 % 1 reduced distribution, 3rd third of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[2 3];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  
  c_eval('ePara = dmpaB?slow.resample(eDist).norm;',ic)
  c_eval('ePerp1 = ePara.cross(dslE?slow.resample(eDist)).norm;',ic)
  ePerp2 = ePara.cross(ePerp1).norm;
  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
  %axis(hca,'equal')
  axis(hca,'square')
end

for ih = 1:numel(h2)
  h2(ih).FontSize = fontsize;
  %h2(ih).XGrid = 'on';
  %h2(ih).YGrid = 'on';
end
%for ih = 4:6, h2(ih).CLim = [-15 -10.2]; end

%% For all the intervals, plot all f's and average
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?.resample(gseVe?);',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
%c_eval('ePitch = ePitch?;',ic)

h2 = setup_subplots(2,1);
%h2().Position(2) = 0.2;
%h2.Position(4) = 0.7;

times_phi_utc = ef1D.tlim(tint_phi).time.utc; times_phi_utc = times_phi_utc(:,12:23);

fontsize = 12;
  
colors = mms_colors('xzy');
isub = 1;
%ndist_to_combine = 3;
n_int = 5;
lowerelim_peaks = 400;
Elim_axis = 4000;
colors = mms_colors('matlab');
colors = colors(3:end,:);
clear plot_f_phi f_sep_utc hmark_all str_Emax Emax Emax_within_interval 
clear Emax_av Emax_av_ind plot_f_phi_Emax Emax_utc f_phi_utc

if 1 % divide f into intervals, common for both v,E x-axis, mark TSeries
  fred_phi = ef1D.tlim(tint_phi);
  plot_v_phi = mean(fred_phi.depend{1},1);
  plot_E_phi = sign(plot_v_phi)*0.5*units.me.*(plot_v_phi*1e3).^2/units.e;   
  
  fred_phi_peaks = fred_phi;   
  fred_phi_peaks.data(abs(repmat(plot_E_phi,fred_phi_peaks.length,1))<lowerelim_peaks) = 0;
  [val,ind] = max(fred_phi_peaks.data');
  Emax = (plot_E_phi(ind));
     
  nf = fred_phi.length;
  [ind1,ind2] = ind12(nf,n_int);
  for i_int = 1:n_int
    f_phi_utc{i_int,1} = sprintf('%s - %s',times_phi_utc(ind1(i_int),:),times_phi_utc(ind2(i_int),:));
    Emax_tmp = Emax(ind1(i_int):ind2(i_int));
    [val_tmp_Emax,ind_tmp_Emax] = max(abs(Emax_tmp));
    Emax_within_interval(i_int) = Emax_tmp(ind_tmp_Emax);
    time_Emax_within_interval{i_int} = fred_phi_peaks.time(ind1(i_int)-1+ind_tmp_Emax);
    %Emax_utc = 
        
    plot_f_phi(i_int,:) = mean(fred_phi.data(ind1(i_int):ind2(i_int),:),1);    
    plot_f_phi_Emax(i_int,:) = fred_phi.data(ind1(i_int)+ind_tmp_Emax-1,:);
    
    Emax_tmp = plot_f_phi(i_int,:); Emax_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(Emax_tmp);
    Emax_av(i_int) = (plot_E_phi(ind));
    Emax_av_ind(i_int) = ind;
    f_at_Emax_av(i_int) = plot_f_phi(i_int,ind);
    
    str_Emax{i_int,1} = sprintf('phimax =[%.0f, %.0f] eV',Emax_av(i_int),Emax_within_interval(i_int));
    
    if 0
    [hmark,cc] = irf_pl_mark(h(npanels_large + 2),fred_phi.time([ind1(i_int) ind2(i_int)]).epochUnix');    
    if isa(hmark,'matlab.graphics.primitive.Line')
      hmark.Color = colors(i_int,:);
    elseif isa(hmark,'matlab.graphics.primitive.Patch')
      hmark.YData = hmark.YData + 120;
      hmark.ZData = [1 1 1 1];
      hmark.EdgeColor = [0 0 0];
      hmark.FaceColor = colors(i_int,:);
    end    
    hmark_all{i_int} = hmark;
    end
  end  
end
if 1 % all fred during sep interval, vs +-Epar, max within interal
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_E_phi,plot_f_phi);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    f_tmp = plot_f_phi_Emax(i_int,:);
    f_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(abs(f_tmp));
    h_fered_star = plot(hca,Emax_av(i_int),f_at_Emax_av(i_int),'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'E_{||} (eV)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  hca.XLim = Elim_axis*[-1 1]+ 0.5*mean(Emax);
  hca.XTick = -30000:1000:30000;
  if hca.XLim(2)<max(Emax)
    hca.XLim = max(Emax)*1*[-1 1] + 0.5*max(Emax);
    hca.XTick = -30000:2000:30000;
  end
     
  hca.Title.String = sprintf('E > %g eV',lowerelim_peaks);
  if 0 % plot lowerelim_peaks
    hold(hca,'on')
    plot(hca,lowerelim_peaks*[-1 -1],hca.YLim,'--k',lowerelim_peaks*[1 1],hca.YLim,'--k')
    hold(hca,'off')
  end
  
  set(hca,'colororder',colors)
  if 0
    if plot_E_phi(ind) > 0
      irf_legend(hca,str_Emax,[0.02 0.98],'fontsize',fontsize)
    else
      irf_legend(hca,str_Emax,[0.98 0.98],'fontsize',fontsize)
    end
  end
  legs_Emax = {};
  for i_int = 1:n_int
    time_utc_tmp = time_Emax_within_interval{i_int}.utc;
    %legs_Emax{i_int,1} = sprintf('%s',time_utc_tmp(12:23));
    legs_Emax{i_int,1} = sprintf('%s: acc. pot. = %4.0f eV',time_utc_tmp(12:23),round(Emax_av(i_int)/100)*100);
  end
  irf_legend(hca,legs_Emax,[0.98 0.98],'fontsize',fontsize)
  %hca.Title.String = {'Max phi in indicated tint',sprintf('phimax with lower limit E > %g eV',lowerelim_peaks)};
  hca.Title.String = '';
  hca.YLim = [0 0.6]*1e-3;
  hca.XLim = [-2000 6000];  
end

if 1 % all fred during sep interval, vs +-Epar, max within interal
  hca = h2(isub); isub = isub + 1;
  
  h_fered = plot(hca,plot_E_phi,plot_f_phi_Emax);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)
  hold(hca,'on')
  for i_int = 1:n_int
    f_tmp = plot_f_phi_Emax(i_int,:);
    f_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(abs(f_tmp));
    h_fered_star = plot(hca,plot_E_phi(ind),val,'*','color',colors(i_int,:));
  end
  hold(hca,'off')
  
  hca.XLabel.String = 'E_{||} (eV)';
  hca.YLabel.String = sprintf('f_e (%s)',fred_phi.units);
  hca.Title.String = '';           
  
  set(hca,'colororder',colors)  
  legs_Emax = {};
  for i_int = 1:n_int
    time_utc_tmp = time_Emax_within_interval{i_int}.utc;
    %legs_Emax{i_int,1} = sprintf('%s',time_utc_tmp(12:23));
    legs_Emax{i_int,1} = sprintf('%s: acc. pot. = %g eV',time_utc_tmp(12:23),round(Emax_within_interval(i_int)/100)*100);
  end
  irf_legend(hca,legs_Emax,[0.98 0.98],'fontsize',fontsize) 
  
  hca.YLim = [0 0.6]*1e-3;
  hca.XLim = [-2000 6000];
  hca.XTick = -30000:1000:30000;
end


n_intervals_sep = 5;
vgmax = 60000;
vg = -vgmax:2000:vgmax;
%lowerelim = 40;
if 0 % 1 reduced distribution, 1st half of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[0 1];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  
  c_eval('ePara = dmpaB?slow.resample(eDist).norm;',ic)
  c_eval('ePerp1 = ePara.cross(dslE?slow.resample(eDist)).norm;',ic)
  ePerp2 = ePara.cross(ePerp1).norm;
  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  %axis(hca,'equal')
  axis(hca,'square')
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
end
if 0 % 1 reduced distribution, 2nd half of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[1 2];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  

  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  %axis(hca,'equal')
  axis(hca,'square')
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
end
if 0 % 1 reduced distribution, 3rd third of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[2 3];
  
  c_eval('eDist = ePDist?.tlim(tint_tmp);',ic)
  
  c_eval('ePara = dmpaB?slow.resample(eDist).norm;',ic)
  c_eval('ePerp1 = ePara.cross(dslE?slow.resample(eDist)).norm;',ic)
  ePerp2 = ePara.cross(ePerp1).norm;
  
  nMC = 500;  
  tic; ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc 
  %tic; ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc
  %tic; ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'scpot',scpot.resample(eDist),'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC); toc

  [h_surf,h_axis,h_all] = ef2D_parperp1.plot_plane(hca,'tint',tint_tmp);
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
  hca.XLim = vg([1 end])*1e-3;
  hca.YLim = vg([1 end])*1e-3;
  tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
  hca.Title.String = {sprintf('%s - %s',tint_tmp_utc(1,:),tint_tmp_utc(end,:)),...
                      sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
  %axis(hca,'equal')
  axis(hca,'square')
end

for ih = 1:numel(h2)
  h2(ih).FontSize = fontsize;
  %h2(ih).XGrid = 'on';
  %h2(ih).YGrid = 'on';
end
%for ih = 4:6, h2(ih).CLim = [-15 -10.2]; end
