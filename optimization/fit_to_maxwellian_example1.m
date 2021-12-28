% Make fit of N Maxwellians to 1D reduced FPI distribution
%% Load and prepare data
% Spacecraft id
ic = 1;
units = irf_units;
T = 0.03; % electron, 0.15 for ions

% Time interval of event
tint_burst = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tint_burst = tint_burst + [+5 -5]; % using the above edges causes problem with new EDI files because they have different versions that adjoining file

% Time interval for figure
tint_dist = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.630Z');
tint_dist = irf_time('2017-07-06T13:54:05.56Z','utc>EpochTT') + 60*T*[-1 1];
%tint_dist = tint_burst;

% Set up database
localuser = datastore('local','user');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS'); 
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   

% Load data
% Necessary data
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint_burst);',ic);
c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint_burst+[20 0]));',ic)
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint_burst);',ic);

% Optional data, to specify starting guess and compare to results
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint_burst,?);',ic)
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint_burst,?);',ic)
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint_burst,?);',ic);

%% Reduce or rebin distribution function into cartesian grid
eint = [000 40000];
vint = [-Inf Inf];
vg = (-70:1:70)*1e3;

c_eval('eDist = ePDist?.tlim(tint_dist);',ic)
           
scpot = scPot1.resample(eDist);
ePara = dmpaB1.resample(eDist(1)).norm;
ePerp1 = ePara.cross(irf.ts_vec_xyz(ePara.time,repmat([1 0 0],ePara.length,1))).norm;
ePerp2 = ePara.cross(ePerp1).norm;
%lmn = [ePara.data; ePerp1.data; ePerp2.data];

%orient = [ePara.data; ePerp1.data; ePerp2.data];
orient = [ePerp1.data; ePerp2.data; ePara.data];
lowerelim = 40;
nMC = 100e0;

% For 3D/2D/1D fit (can further be reduced to 1D or 2D) (quite slow, buggy?)
%ef3D_perp1perp2par = eDist.elim([40 inf]).rebin('cart',{vg,vg,vg},orient); % Rebins skymap into 3D cartesian grid

% For 2D/1D fit (can further be reduced to 1D)
ef2D_parperp1 = eDist.reduce('2D',ePara,ePerp1,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC);
ef2D_parperp2 = eDist.reduce('2D',ePara,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC);
ef2D_perp1perp2 = eDist.reduce('2D',ePerp1,ePerp2,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'base','cart','nMC',nMC);

% For 1D fit
ef1D_par = eDist.reduce('1D',ePara.data,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); % reduced distribution along B
ef1D_perp1 = eDist.reduce('1D',ePerp1.data,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); % reduced distribution along B
ef1D_perp2 = eDist.reduce('1D',ePerp2.data,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'vg',vg,'nMC',nMC); % reduced distribution along B

% Get some values for initial guess of search function
c_eval('n0 = ne?.resample(eDist);' ,ic)
c_eval('v0 = gseVe?*orient''; v0 = v0.resample(eDist);' ,ic)
c_eval('vpar = gseVe?*ePara.data''; vpar = vpar.resample(eDist);' ,ic) % comparison to results
c_eval('T0 = mms.rotate_tensor(gseTe?,''rot'',orient(1,:),orient(2,:),orient(3,:)); T0 = T0.resample(eDist);',ic) % xx component is par

%% Make fit, 1D, with funFitVDF()
vdf = ef1D_par;
nPop = 2;
tic
[fitdata,ts] = funFitVDF(vdf,'nPop',nPop,'plot',0);
toc
%% Plot results from funFitVDF()
h = irf_plot(6);
hca = irf_panel('vdf_obs');
irf_spectrogram(hca,vdf.specrec('velocity'),'log'); colormap(hca,irf_colormap('waterfall')) 
hca.YLabel.String = 'v^{obs} (km/s)';
hca.YLabel.Interpreter = 'tex';

hca = irf_panel('vdf_fit');
irf_spectrogram(hca,ts.f.specrec('velocity'),'log'); colormap(hca,irf_colormap('waterfall'))
hold(hca,'on'); irf_plot(hca,ts.vd); hold(hca,'off'); 
hca.YLabel.String = 'v^{fit} (km/s)';
hca.YLabel.Interpreter = 'tex';

hca = irf_panel('n_fit');
irf_plot(hca,ts.n)
hold(hca,'on'); hnobs = irf_plot(hca,n0); hold(hca,'off'); 
irf_legend(hca,'obs',[0.98 0.98],'color',hnobs.Color)
hca.YLabel.String = 'n (cm^{-3})';
hca.YLabel.Interpreter = 'tex';

hca = irf_panel('vd_fit');
irf_plot(hca,ts.vd)
hold(hca,'on'); hvobs = irf_plot(hca,vpar); hold(hca,'off'); 
irf_legend(hca,'obs',[0.98 0.98],'color',hvobs.Color)
hca.YLabel.String = 'v (km/s)';

hca = irf_panel('T_fit');
irf_plot(hca,ts.T)
hold(hca,'on'); hTobs = irf_plot(hca,T0.xx); hold(hca,'off'); 
irf_legend(hca,'obs',[0.98 0.98],'color',hTobs.Color)
hca.YLabel.String = 'T (eV)';

hca = irf_panel('cf_fit');
irf_plot(hca,ts.cf)

irf_plot_axis_align(h)
hlinks = linkprop([irf_panel('vdf_obs'),irf_panel('vdf_fit')],{'CLim','YLim'});
irf_zoom(h,'x',tint_dist)
