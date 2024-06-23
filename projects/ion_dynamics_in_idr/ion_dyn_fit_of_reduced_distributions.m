%% Load data
ic = 3;
units = irf_units;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854

% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
mms.db_init('local_file_db','/Users/cno062/Data/MMS');
db_info = datastore('mms_db');

 %
% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',1:4);
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);

%c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)

% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',1:4);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint); toc',ic);
%c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',tint); toc',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
%c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);

% Load spacecraft position
disp('Loading spacecraft position...')
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4);

% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);

% Particle moments
% Skymap distributions
disp('Loading skymaps...')
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_nobg = iPDist?; iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic)
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?_counts.data./iPDistErr?.data).^2;',ic)

% Pressure and temperature
disp('Loading pressure and temperature...'); tic
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic); toc

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

% Density
disp('Loading density...'); tic;
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic); toc

% Velocity
disp('Loading bulk velocities...'); tic
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',1:4)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); toc
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic)
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic); toc

c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s
c_eval('gseVexB? = 1e3*cross(1e3*gseVe?,1e-9*gseB?.resample(gseVe?));',1:4) % mV/m

% EDR signatures
c_eval('facPepp? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''pp'');',ic); % Peperp1 = Peperp2
c_eval('facPeqq? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''qq'');',ic); % Peperp1 and Peperp2 are most unequal
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic);
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic);

% Compute Q and Dng from facPepp
c_eval('Q? = (facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2)./(facPepp?.yy.data.^2+2*facPepp?.yy.data.*facPepp?.xx.data);',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);
c_eval('Dng? = sqrt(8*(facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2))./(facPepp?.xx.data+2*facPepp?.yy.data);',ic);
c_eval('Dng? = irf.ts_scalar(ne?.time,Dng?);',ic);

% Compute agyrotropy Aphi from facPeqq
c_eval('agyro? = 2*(facPeqq?.yy-facPeqq?.zz)/(facPeqq?.yy+facPeqq?.zz); agyro? = agyro?.abs',ic);

% Compute temperature ratio An
c_eval('Temprat? = facPepp?.xx/(facPepp?.yy);',ic);

c_eval('gseEVexB? = gseE?.resample(gseVexB?) + gseVexB?''; gseEVexB?.name = ''E+VexB gse'';',1:4)

c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[gseCurlEvexB,gseDivEvexB,gseEvexB,~,~,~] = c_4_j('gseR?brsttime','gseEVexB?');
gseCurlEvexB = gseCurlEvexB*1e-6*1e9; % nT/s

nsm = 5;
[gseCurlEvexB_smo,gseDivEvexB_smo,gseEvexB_smo,~,~,~] = c_4_j(gseR1brsttime,gseR2brsttime,gseR3brsttime,gseR4brsttime,gseEVexB1.smooth(nsm),gseEVexB2.smooth(nsm),gseEVexB3.smooth(nsm),gseEVexB4.smooth(nsm));
gseCurlEvexB_smo = gseCurlEvexB_smo*1e-6*1e9; % nT/s

disp('Done.')

%% Remove one counts of distributions
c_eval('iPDist?_mm_! = iPDist?.movmean(!,''RemoveOneCounts'',iPDist?_counts);',ic,[5 9])

%% Reduce distributions
vint = [-Inf Inf];
if 1 % ion LMN, ielim
  tint_ifred =  irf.tint('2017-07-11T22:32:00.00Z/2017-07-11T22:35:20.00Z');
  ielim = [500 40000];
  c_eval('if1DL?_elim_mm_! = iPDist?_mm_!.elim(ielim).tlim(tint_ifred).reduce(''1D'',L,''vint'',vint);',ic,9)
  c_eval('if1DM?_elim_mm_! = iPDist?_mm_!.elim(ielim).tlim(tint_ifred).reduce(''1D'',M,''vint'',vint);',ic,9)
  c_eval('if1DN?_elim_mm_! = iPDist?_mm_!.elim(ielim).tlim(tint_ifred).reduce(''1D'',N,''vint'',vint);',ic,9)
end


%% Estimating the ion components
tint_fit =  irf.tint('2017-07-11T22:33:20.00Z/2017-07-11T22:34:30.00Z');
tint_fit =  irf.tint('2017-07-11T22:33:30.00Z/2017-07-11T22:34:20.00Z');
vdf = if1DN3_elim_mm_90.tlim(tint_fit);
%vdf = if1DN1_nobg.tlim(tint_fit);
vdf = vdf(1:1:vdf.length);
%vdf = vdf(fix(vdf.length/2));
%vdf = vdf([1 vdf.length]);
nPop = 3; % Three populations give colder and faster counterstreaming beams
X0 = [0.01e6, -1000e3, 500e3,...
      0.01e6,  1000e3, 500e3,...
      0.01e6,  0,      1000e3];

nPop = 2; % Two populations looks ok, but beams are warmer and slower to cover region around v=0
X0 = [0.01e6, -1000e3, 500e3,...
      0.01e6,  1000e3, 500e3];

if 0
nPop = 4; % Three populations give colder and faster counterstreaming beams
X0 = [0.01e6, -1000e3, 500e3,...
      0.01e6,  1000e3, 500e3,...
      0.01e6,  -1000e3,      1000e3,...
      0.01e6,  1000e3,      1000e3];
end
%tic; [fitdata_,ts_] = funFitVDF(vdf(100),'nPop',nPop,'plot',1,'guessprevious',0,'X0',X0,'weight',repmat([0 1 1],1,nPop)); toc;
tic; [fitdata_mm_9,ts_mm_9_w100] = funFitVDF(vdf,'nPop',nPop,'plot',0,'guessprevious',0,'X0',X0,'weight',repmat([1 0 0],1,nPop)); toc;

%% Plot fit
h = irf_plot(7);
hca = irf_panel('vdf_obs');
irf_spectrogram(hca,vdf.specrec('velocity'),'lin')
hca = irf_panel('vdf_fit');
irf_spectrogram(hca,ts.f.specrec('velocity'),'lin')
hold(hca,'on'); irf_plot(hca,ts.vd,':'); hold(hca,'off'); 
hca = irf_panel('n_fit');
irf_plot(hca,ts.n)
hca = irf_panel('vd_fit');
irf_plot(hca,ts.vd)
hca = irf_panel('T_fit');
irf_plot(hca,ts.T)
hca = irf_panel('T*n_fit');
irf_plot(hca,ts.T.*ts.n)
hca = irf_panel('cost function');
irf_plot(hca,ts.cf)

irf_plot_axis_align(h)
irf_zoom(h,'x',vdf.time([1 vdf.length]))
hlinks = linkprop([irf_panel('vdf_obs'),irf_panel('vdf_fit')],{'CLim','YLim'});
colormap(pic_colors('candy4'))








