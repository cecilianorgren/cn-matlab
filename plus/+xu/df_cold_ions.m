ic = 1;

% Set datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/data/mms');
db_info = datastore('mms_db');   
localuser = datastore('local','user');

%% Tint
event = 1;
switch event
  case 2 % 20170622030103
    tint = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:03:03.00Z');    
    tint_df = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:03:03.00Z');    
  case 1
    tint = irf.tint('2017-06-22T02:32:53.00Z/2017-06-22T02:33:23.00Z');
    tint_df = irf.tint('2017-06-22T02:32:53.00Z/2017-06-22T02:33:23.00Z');
  case 3
    tint = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:03:03.00Z');    
    tint_df = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:03:03.00Z');    
end
tint_utc = tint.utc;

eventPath = ['/Users/cecilia/GoogleDrive/Research/Yin/Cold_ions/' tint_utc(1,[1:10 12:13 15:16 18:19]) '/']; % for printing
mkdir(eventPath)

%% Load data
% Particle distributions: electrons and ions
disp('Loading particle distributions...')
c_eval('tic; ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
c_eval('tic; iPDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)

% Setting tint from ePDist1
c_eval('tint = ePDist?([1 ePDist?.length]).time;',ic)
    
% Magnetic field
disp('Loading magnetic field...')
c_eval('dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',1:4);
c_eval('gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('gseB?scm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',1:4);

% Electric field
disp('Loading electric field...')
c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);

% Load spacecraft position
disp('Loading spacecraft position...')

try
  R = mms.get_data('R_gse',tint);
  if ~all([isfield(R,'gseR1') isfield(R,'gseR1') isfield(R,'gseR2') isfield(R,'gseR3') isfield(R,'gseR4')])  
    % not the right data fielsd, try to load from irfu database instead
    db_info = datastore('mms_db');   
    mms.db_init('local_file_db','/data/mms');
    R = mms.get_data('R_gse',tint);
    mms.db_init('local_file_db',db_info.local_file_db_root);
  end
  if size(R.gseR1,2) == 4
    c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
  else
    c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
  end
end
% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

% Pressure and temperature
disp('Loading pressure and temperature...')
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic)

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

% Density
disp('Loading density...')
c_eval('tic; ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?); toc;',ic);
c_eval('tic; ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?); toc;',ic);

% Velocity
disp('Loading bulk velocities...')
c_eval('tic; gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('tic; gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?); toc;',ic)

%% Prepare data

c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)

% Current: use electron density to calculate current
units = irf_units;
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?);',ic);

% Pitchangle distribution
%c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,13); ePitch? = ePitch?.convertto(''s^3/km^6'');',ic)

c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = ''km/s'';',ic) % km/s

c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.units = ''mV/m'';',ic)
c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.units = ''mV/m'';',ic)

c_eval('gseEVexB? = gseE?.resample(gseVexB?.time)+gseVexB?; gseEVexB?.name = ''E+VexB'';',ic)

c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
c_eval('beta?e = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
c_eval('beta?i = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
c_eval('beta? = beta?i + beta?e.resample(beta?i);',ic)

c_eval('Ptot? = PB?.resample(gsePi?)+gsePi?.trace/3;',ic);

% Current, curlometer
if all([~isempty(gseB1) ~isempty(gseB2) ~isempty(gseB3) ~isempty(gseB4)])
  c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
  [Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
  gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
  gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
  gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';
else
  gseJcurl = irf.ts_vec_xyz(gseB1.time,gseB1.data*NaN);
end

%% Fred
eint = [50 3000];
iDist = iPDist1.tlim(tint).elim(eint);
Vi = gseVi1.tlim(iDist.time).resample(iDist);

scpot = scPot1.resample(iDist);
lowerelim = 0 + scpot*0;
vgmax = 2000;
vg = -vgmax:100:vgmax;

tic; if1D_par = iDist.reduce('1D',dmpaB1.resample(iDist).norm,'vg',vg); toc % reduced distribution along B
tic; if1D_x = iDist.reduce('1D',[1 0 0],'vg',vg); toc % reduced distribution along B
tic; if1D_y = iDist.reduce('1D',[0 1 0],'vg',vg); toc % reduced distribution along B
tic; if1D_z = iDist.reduce('1D',[0 0 1],'lowerelim',lowerelim,'vg',vg); toc % reduced distribution along B
%tic; ef1D_sep = eDist.reduce('1D',dmpaB1.resample(eDist).norm,'scpot',scpot,'lowerelim',lowerelim,'vg',vg); toc % reduced distribution along B
%tic; ef1D_sheet = eDist.reduce('1D',dmpaB1.resample(eDist).norm,'scpot',scpot,'lowerelim',lowerelim,'vg',vg); toc % reduced distribution along B
%tic; ef1D_lobe = eDist.reduce('1D',dmpaB1.resample(eDist).norm,'scpot',scpot,'lowerelim',lowerelim,'vg',vg); toc % reduced distribution along B
disp('Done.')

%% Plot overview figure with focus on ions
npanels = 8;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
cmap = colormap('jet');

fontsize = 10;

iisub = 0;
zoomy = [];

if 1 % B
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Pressures
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Pressure');
  set(hca,'ColorOrder',mms_colors('123'))  
  c_eval('irf_plot(hca,{PB?,facPi?.trace/3,Ptot?},''comp'');',ic)
  hca.YLabel.String = {'Pressure','(nP)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'P_{B}';'P_{i}';'P_{i,tot}'},[1.01 0.9],'fontsize',12);    
  irf_zoom(hca,'y')
end
if 1 % ne
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 0 % J 
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % Vi
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % V ExB
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x.tlim(tint),gseVExB?.y.tlim(tint),gseVExB?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % i DEF omni 64
  iisub = iisub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_i','(eV)'};   
  colormap(hca,cmap) 
  hca.YLim(1) = 100;
end
if 1 % i psd vpar
  iisub = iisub + 1;
  hca = irf_panel('fi ref par');  
  irf_spectrogram(hca,if1D_par.specrec('velocity_1D'));  
  hca.YLim = vgmax*[-1 1];
  hca.YLabel.String = {'v_{i,par}','(km/s)'}; 
  %irf_legend(hca,[num2str(fred_to_plot.ancillary.vint(1),'%.0f') '<v_\perp<' num2str(fred_to_plot.ancillary.vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])  
end
if 0 % i psd vpar
  iisub = iisub + 1;
  hca = irf_panel('fi ref x');  
  irf_spectrogram(hca,if1D_x.specrec('velocity_1D'));  
  hca.YLim = vgmax*[-1 1];
  hca.YLabel.String = {'v_{i,x}','(km/s)'}; 
  %irf_legend(hca,[num2str(fred_to_plot.ancillary.vint(1),'%.0f') '<v_\perp<' num2str(fred_to_plot.ancillary.vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])  
end
if 0 % i psd vpar
  iisub = iisub + 1;
  hca = irf_panel('fi ref y');  
  irf_spectrogram(hca,if1D_y.specrec('velocity_1D'));  
  hca.YLim = vgmax*[-1 1];
  hca.YLabel.String = {'v_{i,y}','(km/s)'}; 
  %irf_legend(hca,[num2str(fred_to_plot.ancillary.vint(1),'%.0f') '<v_\perp<' num2str(fred_to_plot.ancillary.vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])  
end
if 0 % i psd vpar
  iisub = iisub + 1;
  hca = irf_panel('fi ref z');  
  irf_spectrogram(hca,if1D_z.specrec('velocity_1D'));  
  hca.YLim = vgmax*[-1 1];
  hca.YLabel.String = {'v_{i,z}','(km/s)'}; 
  %irf_legend(hca,[num2str(fred_to_plot.ancillary.vint(1),'%.0f') '<v_\perp<' num2str(fred_to_plot.ancillary.vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])  
end
if 1 % E
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E perp
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x.tlim(tint),gseE?perp.y.tlim(tint),gseE?perp.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % E par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

irf_zoom(h,'x',tint)

irf_zoom(h(zoomy),'y') 
irf_plot_axis_align

h(1).Title.String = sprintf('MMS %g',ic);

for ii = 1:npanels
  h(ii).FontSize = 12;
end

for ii = 7:10,  h(ii).YLim = 900*[-1 1]; %h(ii).CLim = [-2 -0];
end

%% Plot electron particle distributions
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?.resample(gseVe?);',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
%c_eval('ePitch = ePitch?;',ic)

times_phi_utc = ef1D.tlim(tint_phi).time.utc; times_phi_utc = times_phi_utc(:,12:23);

fontsize = 8;
  
colors = mms_colors('xzy');
isub = 1;
if 1 % average fred during sep, lobe and sheet interval
  hca = h2(isub); isub = isub + 1;
  
  fered_sep = ef1D.tlim(tint_sep);
  plot_v_sep = nanmean(fered_sep.depend{1},1);
  plot_f_sep = nanmean(fered_sep.data,1);
  
  fered_lobe = ef1D.tlim(tint_lobe);
  plot_v_lobe = nanmean(fered_lobe.depend{1},1);
  plot_f_lobe = nanmean(fered_lobe.data,1);
  
  fered_sheet = ef1D.tlim(tint_sheet);
  plot_v_sheet = nanmean(fered_sheet.depend{1},1);
  plot_f_sheet = nanmean(fered_sheet.data,1);
  
  h_fered = plot(hca,plot_v_lobe*1e-3,plot_f_lobe,...
                     plot_v_sheet*1e-3,plot_f_sheet,...
                     plot_v_sep*1e-3,plot_f_sep);
  c_eval('h_fered(?).Color = colors(?,:);',1:3)
  hca.XLabel.String = 'v (10^3 km/s)';
  hca.YLabel.String = sprintf('f_e (%s)',fered_sep.units);
  hca.XLim = [-40 40];
  
  set(hca,'colororder',colors);
  tint_str = {sprintf('%s - %s',tint_lobe_utc(1,:),tint_lobe_utc(2,:));...
              sprintf('%s - %s',tint_sheet_utc(1,:),tint_sheet_utc(2,:));...
              sprintf('%s - %s',tint_sep_utc(1,:),tint_sep_utc(2,:))};
  irf_legend(hca,tint_str,[0.01 0.98],'fontsize',fontsize)
  hca.Title.String = sprintf('E > %g eV',lowerelim);
end

%ndist_to_combine = 3;
n_int = 3;
lowerelim_peaks = 400;
Elim_axis = 4000;
colors = mms_colors('matlab');
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
    %Emax_utc = 
        
    plot_f_phi(i_int,:) = mean(fred_phi.data(ind1(i_int):ind2(i_int),:),1);    
    plot_f_phi_Emax(i_int,:) = fred_phi.data(ind1(i_int)+ind_tmp_Emax-1,:);
    
    Emax_tmp = plot_f_phi(i_int,:); Emax_tmp(abs(plot_E_phi)<lowerelim_peaks) = 0;
    [val,ind] = max(Emax_tmp);
    Emax_av(i_int) = (plot_E_phi(ind));
    Emax_av_ind(i_int) = ind;
    f_at_Emax_av(i_int) = plot_f_phi(i_int,ind);
    
    str_Emax{i_int,1} = sprintf('phimax =[%.0f, %.0f] eV',Emax_av(i_int),Emax_within_interval(i_int));
    
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
if 1 % all fred during sep interval, vs vpar
  hca = h2(isub); isub = isub + 1;
  h_fered = plot(hca,plot_v_sep*1e-3,plot_f_phi);
  c_eval('h_fered(?).Color = colors(?,:);',1:n_int)  
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
  if plot_E_phi(ind) > 0
    irf_legend(hca,str_Emax,[0.02 0.98],'fontsize',fontsize)
  else
    irf_legend(hca,str_Emax,[0.98 0.98],'fontsize',fontsize)
  end
  hca.Title.String = {'Max phi in indicated tint',sprintf('phimax with lower limit E > %g eV',lowerelim_peaks)};
end

n_intervals_sep = 2;
vgmax = 60000;
vg = -vgmax:2000:vgmax;
%lowerelim = 40;
if 1 % 1 reduced distribution, 1st half of tint_sep
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
if 1 % 1 reduced distribution, 2nd half of tint_sep
  hca = h2(isub); isub = isub + 1;
  tint_tmp = tint_phi(1) + (tint_phi(2)-tint_phi(1))/n_intervals_sep*[1 2];
  
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
  