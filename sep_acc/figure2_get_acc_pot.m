units = irf_units;
ic = 1;
%% Time intervals
% burst interval
tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');

% acceleration channels
events = [3 4 5];
nevents = numel(events);
for ievent = 1:3
  event = events(ievent);
  sep.get_tints;
  tints_phi{ievent} = tint_phi;
  tints_sep{ievent} = tint_sep;
  tints_lobe{ievent} = tint_lobe;
  tints_sheet{ievent} = tint_sheet;
end

% 3 acceleration channels
tint_fred = irf.tint('2017-07-06T00:54:10.00Z/2017-07-06T00:54:31.00Z'); 

%% Load data
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)
c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
%c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
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

c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);

c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)

c_eval('vte?par = irf.ts_scalar(facTe?.time,sqrt(2*units.e*facTe?.xx.data/units.me)*1e-3);',ic)

c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)

c_eval('beta?e = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
c_eval('beta?i = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
c_eval('beta? = beta?i + beta?e.resample(beta?i);',ic)

%% Pick out parameters in the lobe, sheet, and separatrix intervals
c_eval('B_lobe{!} = mean(gseB?.abs.tlim(tints_lobe{!}).data,1);',ic,1:nevents)

c_eval('n_lobe{!} = mean(ne?.tlim(tints_lobe{!}).data,1);',ic,1:nevents)
c_eval('n_sheet{!} = mean(ne?.tlim(tints_sheet{!}).data,1);',ic,1:nevents)
c_eval('n_sep{!} = mean(ne?.tlim(tints_sep{!}).data,1);',ic,1:nevents)

c_eval('be_lobe{!} = mean(beta?e.tlim(tints_lobe{!}).data,1);',ic,1:nevents)
c_eval('be_sheet{!} = mean(beta?e.tlim(tints_sheet{!}).data,1);',ic,1:nevents)
c_eval('be_sep{!} = mean(beta?e.tlim(tints_sep{!}).data,1);',ic,1:nevents)

c_eval('Tepar_lobe{!} = mean(facTe?.tlim(tints_lobe{!}).data(:,1,1),1);',ic,1:nevents)
c_eval('Tepar_sheet{!} = mean(facTe?.tlim(tints_sheet{!}).data(:,1,1),1);',ic,1:nevents)
c_eval('Tepar_sep{!} = mean(facTe?.tlim(tints_sep{!}).data(:,1,1),1);',ic,1:nevents)

c_eval('Teperp_lobe{!}  = 0.5*mean(facTe?.tlim(tints_lobe{!}).data(:,2,2)+facTe?.tlim(tints_lobe{!}).data(:,3,3),1);',ic,1:nevents)
c_eval('Teperp_sheet{!} = 0.5*mean(facTe?.tlim(tints_sheet{!}).data(:,2,2)+facTe?.tlim(tints_sheet{!}).data(:,3,3),1);',ic,1:nevents)
c_eval('Teperp_sep{!}   = 0.5*mean(facTe?.tlim(tints_sep{!}).data(:,2,2)+facTe?.tlim(tints_sep{!}).data(:,3,3),1);',ic,1:nevents)

c_eval('Te_lobe{?} = (2*Teperp_lobe{?}+Tepar_lobe{?})/3',1:nevents);
if 0
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
end
 
%% Reduced electron distribution
ic = 1;
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

ef1D1_orig = ef1D_orig;
ef1D1_nobg = ef1D_nobg;

if 1
  ef1D1_nobg_smooth3 = ef1D1_nobg.smooth(3);
  ef1D1_nobg_smooth6 = ef1D1_nobg.smooth(6);
else
  ef1D1_nobg1 = ef1D1_nobg.resample(ef1D1_nobg.time(1:3:end));
  ef1D1_nobg2 = ef1D1_nobg.resample(ef1D1_nobg.time(2:3:end));
  ef1D1_nobg3 = ef1D1_nobg.resample(ef1D1_nobg.time(3:3:end));
  ef1D1_nobg.data(1:3:end,:) = ef1D1_nobg1.data;
  ef1D1_nobg.data(2:3:end,:) = ef1D1_nobg2.data;
  ef1D1_nobg.data(3:3:end,:) = ef1D1_nobg3.data;
end
%ef1D_orig = ef1D;
%ef1D = ef1D_new;
if 0
  %%
  fsheet_orig = ef1D_orig.tlim(tint_sheet);
  fsheet_nobg = ef1D_nobg.tlim(tint_sheet);
  plot(fsheet_orig.depend{1}(1,:),mean(fsheet_orig.data,1),...
       fsheet_nobg.depend{1}(1,:),mean(fsheet_nobg.data,1))
end

%% Get acceleration potential
f_lobe = ef1D1_nobg.tlim(tint_fred(1) + [0 1]); 
f_lobe = max(mean(f_lobe.data,1));  
f_prom = 0.1*f_lobe*n_sep{1}/n_lobe{1};

eint = [000 40000];
quadrant = [0 0; 1 0];
if any(quadrant([1 4])) % top left or bottom right, inflow is antiparallel to B
  vint = [-100000 0];
else  % top right or bottom left, inflow is parallel to B
  vint = [0 100000];
end
relativeminpeakprominence = 0.10;
minpeakprominence = f_prom;
resample_timestep = 3;

doPlot = 1;
% Need to make for-loop for printing, since figure is done within script
for iic = ic          
  figure(11)
  c_eval('[tsAccPot?_nobg_rel,tmppPeakF,tmpTsFe] = find_acc_pot(ef1D?_nobg.smooth(3),''eint'',eint,''vint'',vint,''relativeminpeakprominence'',relativeminpeakprominence,''resample'',resample_timestep,''plot'',doPlot);',iic)
  c_eval('tsAccVel?_nobg_rel = irf.ts_scalar(tsAccPot?_nobg_rel.time,sqrt(tsAccPot?_nobg_rel.data*units.e*2/units.me));',iic)
  hcf = gcf;   
  
  hcf.Position = [0 100 700 800];
  %cn.print(sprintf('event%g_tsAccPot%g_nobg_rel',event,iic),'path',printAccPotPath)

  figure(12)
  c_eval('tsAccPot?_nobg_abs = find_acc_pot(ef1D?_nobg,''eint'',eint,''vint'',vint,''minpeakprominence'',minpeakprominence,''resample'',resample_timestep,''plot'',doPlot);',iic)    
  c_eval('tsAccVel?_nobg_abs = irf.ts_scalar(tsAccPot?_nobg_abs.time,sqrt(tsAccPot?_nobg_abs.data*units.e*2/units.me));',iic)
  hcf = gcf;    
  hcf.Position = [700 100 700 800];
  %cn.print(sprintf('event%g_tsAccPot%g_nobg_abs',event,iic),'path',printAccPotPath)

  c_eval('[tmpTsAccPot?,tmppPeakF?,tmpTsFe?] = find_acc_pot(cond.dist,''eint'',eint,''vint'',vint,''minpeakheight'',minpeakheight,''minpeakprominence'',minpeakprominence,''resample'',3);',ic)
  
end
 
%% Plot
hcf = figure(13);
%hcf.Position = [1400 100 700 800];
h = irf_plot(1);
h(1).Position(1) = [0.14];
h(1).Position(3) = [0.7];

ylim = [-50 50];
isub = 1;
if 1 % original fred 
  hca = irf_panel('fred');
  fred_min = 1e-6;
  fred_to_plot = ef1D1_nobg; 
  fred_to_plot.data(fred_to_plot.data < fred_min) = NaN;
  [~,hcb] = irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));    

  
  if 1 % plot acc. vel
    hold(hca,'on')    
    plot_vel = tsAccVel1_nobg_abs*1e-6;
    for itint = 1:3
      hleg(1) = irf_plot(hca,plot_vel.tlim(tints_sep{itint}),'-','linewidth',3,'markersize',15,'color',0*[0.5 0.5 0.5].^1.5);
    end    
    hold(hca,'off')
  end
  if 1 % plot el. vel
    hold(hca,'on')  
    plot_v = gseVe1par*1e-3;    
    hleg(2) = irf_plot(hca,plot_v.resample(plot_v.time(1:3:end)),'k-','linewidth',1);        
    hold(hca,'off')
  end
  if 1 % plot thermal velocity interval
    hold(hca,'on')        
    plot_v = gseVe1par*1e-3-vte1par*1e-3;
    hleg(3) = irf_plot(hca,plot_v.resample(plot_v.time(1:3:end)),'k--','linewidth',1);
    plot_v = gseVe1par*1e-3+vte1par*1e-3;
    irf_plot(hca,plot_v.resample(plot_v.time(1:3:end)),'k--','linewidth',1);
    hold(hca,'off')
  end
  for itint = []%1:3
    irf_pl_mark(hca,tints_phi{itint},'k')
  end
  if 1 % mark potential in text
    % interval 1   
    clear phi vel
    for itint = 1:3            
      ind = find(tsAccPot1_nobg_abs.tlim(tints_sep{itint}).data*1e-6 == max(tsAccPot1_nobg_abs.tlim(tints_sep{itint}).data*1e-6));
      ind = ind(1);      
      phi{itint} = tsAccPot1_nobg_abs.tlim(tints_sep{itint}).data(ind);
      vel{itint} = tsAccVel1_nobg_abs.tlim(tints_sep{itint}).data(ind)*1e-6;
    end
    fontsize = 14;
    irf_legend(hca,{sprintf('%.0f V',round(phi{1}/100)*100)},[0.2 1.02],'fontsize',fontsize,'color',[0 0 0])
    irf_legend(hca,{sprintf('%.0f V',round(phi{2}/100)*100)},[0.63 1.02],'fontsize',fontsize,'color',[0 0 0])
    irf_legend(hca,{sprintf('%.0f V',round(phi{3}/100)*100)},[0.86 1.02],'fontsize',fontsize,'color',[0 0 0])
    %tx = tsAccPot1_nobg_abs.tlim(tints_sep{itint}).time(ind);
    %xx = [tx tx]-tint(1);
    %yy = ([ylim(2) vel]-ylim(1))/diff(ylim);
    %annotation('textarrow',xx/diff(hca.XLim),yy,'String',sprintf('%.0f eV',round(phi/100)*100))

  end
  if 0 % plot vlim
    hold(hca,'on')
    irf_plot(hca,irf.ts_scalar(fred.time([1 end]),vlim(1)*[1 1])*1e-3,'--k');
    irf_plot(hca,irf.ts_scalar(fred.time([1 end]),vlim(2)*[1 1])*1e-3,'--k');
    hold(hca,'off')
  end
  
  hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
  if 0%plotAccVel
    hlegs = legend(h_,legends,'location','best');
    hlegs.Title.String = 'downsampled to:';
  end
  legend(hleg,{'v_{acc}','v_{e||} (bulk)','v_{e,||}\pm v_{te,||}'},'box','on','location','southeast')
  colormap(hca,'jet')
  hca.YLim = [-50 50];
  hca.FontSize = 14;
  hca.Position(1) = 0.12;
  hca.Position(3) = 0.75;
  hca.Position(2) = 0.2;
  hca.Position(4) = 0.7;
  
  hcb.Position(1) = 0.89;
  hcb.Position(2) = 0.2;
  hcb.Position(4) = 0.7;
  hca.YTick = [-60:20:60];
end
