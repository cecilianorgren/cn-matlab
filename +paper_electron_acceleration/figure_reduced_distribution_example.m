ic = 1:4;

%% Set tint
tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');
tint_fig = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:54:45.00Z');
tint_lobe = irf.tint('2017-07-06T00:54:07.00Z/2017-07-06T00:54:08.00Z');
tint_sheet = irf.tint('2017-07-06T00:54:18.70Z/2017-07-06T00:54:19.50Z');
tint_sep = irf.tint('2017-07-06T00:54:13.50Z/2017-07-06T00:54:16.50Z');
tint_phi = irf.tint('2017-07-06T00:54:14.00Z/2017-07-06T00:54:15.50Z');
tint_fred = tint_phi;
times_fred = irf_time('2017-07-06T00:54:14.20Z','utc>epochtt') + [0 0.1 0.2 0.3];
  
%% Load datastore
mms.db_init('local_file_db','/Volumes/Nexus/data');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
pathLocalUser = ['/Users/' localuser '/'];

%% Load data
units = irf_units;

% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);

% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('tic; scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);

% Skymap distributions
disp('Loading skymaps...')
%c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
%c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('tic; [iPDist?,iPDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)
c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)

% Remove all one-count "noise"
%c_eval('iPDist?.data(iPDist?.data<iPDistErr?.data*1.1) = 0;',ic)
%c_eval('ePDist?.data(ePDist?.data<ePDistErr?.data*1.1) = 0;',ic)

disp('Done loading data.');

%% Make reduced distribution, remove background
ZeroNaN = 0;
tic; [eDist1_bgremoved, eDist1_bg, ephoto_scale] = mms.remove_edist_background(ePDist1, 'tint', tint_fred, 'ZeroNaN', ZeroNaN); toc
eDist = eDist1_bgremoved;  
eDist_bgremoved = eDist1_bgremoved;  
eDist_original = ePDist1;

%% Make reduced distribution, reduce distributon
% scpot_lim = scPot1.resample(eDist);
% lowerelim = 0;
% lowerelim_ts = lowerelim  + 0*scpot_lim;
% eLine = dmpaB1.resample(eDist).norm;
% energies = sqrt(2*eDist.depend{1}(1,:)*units.eV/units.me)/1000; % km/s
% vg = [-energies(end:-1:1) 0 energies];
% vg(abs(vg)>80000) = [];
% tic; ef1D = eDist.reduce('1D',eLine,'vint',vint,'scpot',scpot_lim,'lowerelim',lowerelim_ts,'vg',vg); toc % reduced distribution along B

%% Plot acceleration potential and electron particle distributions
nrows = times_fred.length;
ncols = 3;
npanels = nrows*ncols;
for ipanel = 1:npanels
  h2(ipanel) = subplot(nrows,ncols,ipanel);  
end
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?.resample(gseVe?);',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
  
colors = mms_colors('xzy');
isub = 1;

for it = 1:times_fred.length
  tint_tmp = times_fred(it) + 3*0.99*0.5*0.030*[-1 1];
  edist_bgremoved_tmp = eDist_bgremoved.tlim(tint_tmp);
  edist_original_tmp = eDist_original.tlim(tint_tmp);
  
  ePara = dmpaB1slow.resample(edist_original_tmp).norm;
  ePerp1 = ePara.cross(dslE1slow.resample(edist_original_tmp)).norm;
  ePerp2 = ePara.cross(ePerp1).norm;

  nMC = 700;
  scpot = scPot1.resample(eDist).tlim(tint_tmp);
  lowerelim = 000;
  lowerelim_ts = lowerelim  + 0*scpot;
  energies = sqrt(2*edist_original_tmp.depend{1}(1,:)*units.eV/units.me)/1000; % km/s
  vg = [-energies(end:-1:1) 0 energies];
  vg(abs(vg)>80000) = [];
  vg2d = -80000:2000:80000;
  
  e1d_original = edist_original_tmp.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim_ts,'vg',vg,'nMC',nMC); % reduced distribution along B
  e1d_bgremoved = edist_bgremoved_tmp.reduce('1D',ePara,'vint',vint,'scpot',scpot,'lowerelim',lowerelim_ts,'vg',vg,'nMC',nMC); % reduced distribution along B
  
  e2d_original = edist_original_tmp.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim_ts,'vg',vg2d,'base','cart','nMC',nMC); % reduced distribution along B
  e2d_bgremoved = edist_bgremoved_tmp.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim_ts,'vg',vg2d,'base','cart','nMC',nMC); % reduced distribution along B
  
  if 1 % 1d, original and bg removed
    hca = h2(isub); isub = isub + 1;        
    h_fered = plot(hca,mean(sign(e1d_original.depend{1})*0.5*units.me.*(e1d_original.depend{1}*1e3).^2/units.e,1),[nanmean(e1d_original.data,1); nanmean(e1d_bgremoved.data,1)]);
    hca.XLabel.String = 'E (eV)';
    hca.YLabel.String = sprintf('f_e (%s)',ef1D_tmp.units);
    hca.XLim = 2000*[-1 1];    
    irf_legend(hca,{sprintf('E>%g eV',lowerelim);sprintf('scpot = %.0f V',nanmean(scpot.data))},[0.02 0.98])    
    hca.Title.String = sprintf('%s, ndist = %g, ZeroNaN = %g',irf_time(times_fred(it),'epochtt>utc_yyyy-dd-mmTHH:MM:SS.mmm'),e1d_original.length,ZeroNaN);
  end
  if 1 % 2d reduced distribution, original
    hca = h2(isub); isub = isub + 1;    
            
    [h_surf,h_axis,h_all] = e2d_original.plot_plane(hca);
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
    hca.XLim = vg([1 end])*1e-3;
    hca.YLim = vg([1 end])*1e-3;
    tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
    hca.Title.String = {sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
    %axis(hca,'equal')
    axis(hca,'square')
    hca.Title.FontSize = fontsize;
    hca.FontSize = fontsize;
  end
  if 1 % 2d reduced distribution, bg removed
    hca = h2(isub); isub = isub + 1;    
            
    [h_surf,h_axis,h_all] = e2d_bgremoved.plot_plane(hca);
    hca.XLabel.String = 'v_{||} (10^3 km/s)';
    hca.YLabel.String = 'v_{\perp1} (10^3 km/s)';    
    hca.XLim = vg([1 end])*1e-3;
    hca.YLim = vg([1 end])*1e-3;
    tint_tmp_utc = eDist.time.utc; tint_tmp_utc = tint_tmp_utc(:,12:23);
    hca.Title.String = {sprintf('<v_{perp,1}> = [%.2f,%.2f,%.2f]',mean(ePerp1.data,1))};
    %axis(hca,'equal')
    axis(hca,'square')
    hca.Title.FontSize = fontsize;
    hca.FontSize = fontsize;
  end
end
colormap('jet')

%% Make ts plot
ic = 1;
figure(204)
npanels = 5;
[h,h2] = initialize_combined_plot(npanels,3,2,0.4,'vertical');
iisub = 0;
cmap = colormap('jet');
zoomy = [];
fred_min = 1e-6;
fontsize = 12;

if 1 % e DEF omni 64
  iisub = iisub + 1;
  hca = irf_panel('eDEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
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
if 1 % e psd vpar
  iisub = iisub + 1;
  hca = irf_panel('fe reduced');
  fred_to_plot = ef1D; fred_to_plot.data(fred_to_plot.data < fred_min) = NaN;
  irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));  
  hca.YLim = fred_to_plot.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  %irf_legend(hca,[num2str(fred_to_plot.ancillary.vint(1),'%.0f') '<v_\perp<' num2str(fred_to_plot.ancillary.vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  %irf_legend(hca,['E_{e} >' num2str(lowerelim) ' eV'],[0.98 0.99],'color',0*[1 1 1],'fontsize',fontsize)
  %irf_legend(hca,['f_{e} >' num2str(fred_min) ' s/m^4'],[0.98 0.70],'color',0*[1 1 1],'fontsize',fontsize)
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
if 1 % E par
  iisub = iisub + 1;
  zoomy = [zoomy iisub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end

irf_zoom(h,'x',tint_fred)
irf_zoom(h(zoomy),'y') 
irf_plot_axis_align
h(1).Title.String = sprintf('MMS %g',ic);
hca = irf_panel('eDEF omni');  hca.XGrid = 'off'; hca.YGrid = 'off';
