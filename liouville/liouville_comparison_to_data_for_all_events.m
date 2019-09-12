%% Load data, density and bulk velocity map
[PSI,VPSI,map_N,map_V] = get_map_liouville_moments; % n and v

%% Specify event
event = 1;

%% Load specfic data defined for that event 
% Could be potential shape, eye fit parameters, etc.

%% Load data, specific to event
doLoadData = 0;
if doLoadData
  ic = 1;
  tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');
  units = irf_units;

  sep.get_tints;

  mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
  db_info = datastore('mms_db');   
  localuser = datastore('local','user');
  pathLocalUser = ['/Users/' localuser '/'];

  c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
  c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);

  c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
  c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
  c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)

  c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
  c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)

  tint_fred = irf.tint('2017-07-06T08:16:37.00Z/2017-07-06T08:16:41.00Z');
  c_eval('eDist = ePDist?.tlim(tint_fred);',ic)

  % Remove background
  nSecondary = 5;
  nPhoto = 1;
  tic; eDist_nobg = mms.remove_edist_background(eDist,'nSecondary',nSecondary,'Nphotoe_art',nPhoto); toc;

  % Reduced electron distribution
  eint = [00 40000];
  lowerelim = 000;
  vint = [-Inf Inf];
  %tint_fred = tint_fred;%tint_phi;


  eDist_orig = eDist;

  c_eval('scpot = scPot?.resample(eDist);',ic)
  c_eval('dir_red = dmpaB?.resample(eDist).norm;',ic)
  energies = sqrt(2*eDist.depend{1}(1,:)*units.eV/units.me)/1000; % km/s
  vgmax = 70000;
  vg = -vgmax:1000:vgmax;
  vg(abs(vg)>70000) = [];

  nMC = 200;
  tic; ef1D_orig = eDist_orig.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'nMC',nMC,'vg',vg); toc % reduced distribution along B
  tic; ef1D_nobg = eDist_nobg.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'nMC',nMC,'vg',vg); toc % reduced distribution along B
end
