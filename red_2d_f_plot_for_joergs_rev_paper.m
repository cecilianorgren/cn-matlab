ic = 1;
units = irf_units;
event = 2;
switch event 
  case 1
    tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');
    tint_fred = irf.tint('2017-07-06T00:54:14.85Z/2017-07-06T00:54:15.00Z');
    %tint_fred = irf.tint('2017-07-06T00:54:22.50Z/2017-07-06T00:54:22.80Z');
    %tint_fred = irf.tint('2017-07-06T00:54:27.20Z/2017-07-06T00:54:27.40Z');
    tint_fred = irf.tint('2017-07-06T00:55:33.70Z/2017-07-06T00:55:33.90Z');
    tint_fred = irf.tint('2017-07-06T00:55:33.90Z/2017-07-06T00:55:34.00Z');
    tint_fred = irf.tint('2017-07-06T00:55:33.60Z/2017-07-06T00:55:33.70Z');
  case 2
    tint = irf.tint('2017-07-03T05:26:13.00Z/2017-07-03T05:27:02.00Z');
    tint_fred = irf.tint('2017-07-03T05:26:49.70Z/2017-07-03T05:26:50.30Z');
    %tint_fred = irf.tint('2017-07-06T00:54:22.50Z/2017-07-06T00:54:22.80Z');
    %tint_fred = irf.tint('2017-07-06T00:54:27.20Z/2017-07-06T00:54:27.40Z'); 
end


mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
eventPath = ['Users/' localuser '/GoogleDrive/PapersInProgress/Paper_ReviewChapter_Acceleration_in_inflow_region_and_along_separatrices/Figures/'];

%% Load data
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);

% Remove all one-count "noise"
c_eval('ePDist?_no1c = ePDist?;',ic)
c_eval('ePDist?_no1c.data(ePDist?_no1c.data<ePDistErr?.data*1.1) = 0;',ic)

%% Remove secondary electron background
ZeroNaN = 0;
Nphoto_art = 0;
nSecondary = 8;
tic; [eDist1_bgremoved, eDist1_bg, ephoto_scale] = mms.remove_edist_background(ePDist1, 'tint', tint_fred, 'ZeroNaN', ZeroNaN, 'Nphotoe_art', Nphoto_art, 'nSecondary', nSecondary); toc

%c_eval('ePDist?_no1c = ePDist?;',ic)
%c_eval('ePDist?_no1c.data(ePDist?_no1c.data<ePDistErr?.data*1.1) = 0;',ic)

c_eval('eDist_long = ePDist?_no1c.tlim(tint_fred+[-2 2]);',ic)
c_eval('eDist_original = ePDist?.tlim(tint_fred);',ic)
c_eval('eDist_bgremoved = eDist?_bgremoved.tlim(tint_fred);',ic)
c_eval('eDist_no1c = ePDist?_no1c.tlim(tint_fred);',ic)

%% Make reduced distribution
c_eval('scpot = scPot?.resample(eDist_original);',ic)
c_eval('scpot_long = scPot?.resample(ePDist?);',ic)
c_eval('par = dmpaB?.resample(eDist_original).norm;',ic)
c_eval('perp1 = dslE?.resample(eDist_original).norm.cross(par); perp1 = perp1.norm;',ic)
c_eval('perp2 = par.cross(perp1); perp2 = perp2.norm;',ic)

vg = -70000:1000:70000;
tic; ef1D_long = eDist_long.reduce('1D',par,'vg',vg,'scpot',scpot_long); toc
tic; ef2D_orig = eDist_original.reduce('2D',par,perp1,'base','cart','vg',vg,'scpot',scpot,'lowerelim',1.5*scpot); toc
tic; ef2D_nobg = eDist_bgremoved.reduce('2D',par,perp1,'base','cart','vg',vg,'scpot',scpot); toc
tic; ef2D_no1c = eDist_no1c.reduce('2D',par,perp1,'base','cart','vg',vg,'scpot',scpot); toc

%% Plot
nrows = 4;
ncols = 1;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols);
isub = 1;

hca = h(isub); isub = isub + 1;
ef1D_long.data(ef1D_long.data<1e-6) = NaN;
irf_spectrogram(hca,ef1D_long.specrec);
irf_pl_mark(hca,tint_fred)

hca = h(isub); isub = isub + 1;
ef2D_orig.plot_plane(hca);
hca.Title.String = 'original';

hca = h(isub); isub = isub + 1;
ef2D_nobg.plot_plane(hca);
hca.Title.String = {'background removed:';sprintf(' Nphoto_art = %g cc, nSecondary = %g cc',Nphoto_art,nSecondary)};
hca.Title.Interpreter ='none';
Nphoto_art = 0;
nSecondary = 2;


hca = h(isub); isub = isub + 1;
ef2D_no1c.plot_plane(hca);
hca.Title.String = 'one-counts removed';


for ipanel = 2:npanels
  axis(h(ipanel),'equal')
  axis(h(ipanel),'square')
  h(ipanel).XLim = 2*0.5*[-70 70];
  h(ipanel).YLim = 2*0.5*[-70 70];
  h(ipanel).XTick = -100:20:100;
  h(ipanel).YTick = -100:20:100;
  h(ipanel).YLabel.String = 'v_{ExB} (10^3 km/s)';
  h(ipanel).XLabel.String = 'v_{||} (10^3 km/s)';
end
colormap('jet')
hlink = linkprop(h(2:end),{'CLim','XLim','YLim'});

irf_timeaxis(h(1));
