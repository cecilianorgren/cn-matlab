%% Load data
% Load datastore
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   
localuser = datastore('local','user');

% Specify time intervals for waves
tint_waves{1} = irf.tint('2017-05-05T20:09:47.900Z/2017-05-05T20:09:47.960Z');
tint_waves{2} = irf.tint('2017-05-05T20:09:49.736Z/2017-05-05T20:09:49.746Z');
tint_waves{3} = irf.tint('2017-05-05T20:09:50.820Z/2017-05-05T20:09:50.839Z');

% Specify time interval for file
ic = 1:4;
tint_search = irf.tint('2017-05-05T20:09:20.00Z/2017-05-05T20:10:00.00Z');
files = mms.db_list_files('mms1_edp_brst_l2_scpot',tint_search);
tint = [files.start files.stop];
%%
listDefatt = mms.db_list_files(['mms', '1', '_ancillary_defatt']);
defatt1 = mms.db_get_variable('mms1_ancillary_defatt','zra',tint);
dataTmp = mms_load_ancillary([listDefatt(end).path, filesep, ...
      listDefatt(end).name], 'defatt');
%zphase1 = irf.ts_scalar()
%%
% Magnetic field
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);

% Electric field
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);

% Spacecraft potential
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);
c_eval('tic; psp?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_psp_brst_l2'',tint); toc;',ic);

% Load scpot dataobj
files_scpot = mms.db_list_files('mms1_edp_brst_l2_scpot',tint);
dobj_scpot = dataobj([files_scpot.path filesep files_scpot.name]);

% Load defatt
c_eval('zphase? = mms.db_get_variable(''mms?_ancillary_defatt'',''zphase'',tint);',ic);
c_eval('zphase? = irf.ts_scalar(irf_time(zphase.time,''ttns>EpochTT''),zphase.zphase);',ic);

%c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic);
%c_eval('defatt? = mms_removerepeatpnts(defatt?);',ic)


% Specify paths for figures and matlab files
fileName = dmpaB1.userData.GlobalAttributes.Logical_file_id;
fileNameSplit = strsplit(fileName{1},'_'); 
eventPath = ['/Users/' localuser '/Research/Events/Hwang2021_probe_timing/'];  
matlabPath = ['/Users/' localuser '/MATLAB/cn-matlab/mms' fileNameSplit{5} '/'];

%% Plot probe data
ic = 1;
npanels = 4;
h = irf_plot(npanels);

hca = irf_panel('E');
c_eval('irf_plot(hca,dslE?);',ic)
irf_legend(hca,{'x','y','z'},[0.98 0.98]);
hca.YLabel.String = {'E DSL','(mV/m)'};

hca = irf_panel('scpot');
c_eval('irf_plot(hca,scPot?);',ic)
hca.YLabel.String = {'V_{sc}','(V)'};

hca = irf_panel('psp');
c_eval('irf_plot(hca,psp?);',ic)
irf_legend(hca,{'1','2','3','4','5','6'},[0.98 0.98]);
hca.YLabel.String = {'Probe2Sc pot','(V)'};
irf_legend(hca,'Probe to spacecraft potential, averaged',[0.02 0.02],'color',[0 0 0])

hca = irf_panel('dcv');
c_eval('irf_plot(hca,dcv?);',ic)
irf_legend(hca,{'1','2','3','4','5','6'},[0.98 0.98]);
hca.YLabel.String = {'Probe2Sc pot','(V)'};
irf_legend(hca,'Individual probes. P1=V1, P2=(V1-0.120*E12), P3=V3, P4=(V3-0.120*E34), P5=V5, P6=(V5-0.0292*E56)',[0.02 0.02],'color',[0 0 0])

c_eval('irf_pl_mark(h(?),tint_waves{!},[0 0 0]);',1:npanels,1:numel(tint_waves))

tint_zoom = [tint_waves{1}(1) tint_waves{end}(2)] + [-2 2];
irf_zoom(h,'x',tint_zoom)

%% Check probe align times for spin plane probes
%mms.probe_align_times
%mms.sc_orient
ic = 1;
c_eval('Exyz = dslE?;',ic)
c_eval('Bxyz = dmpaB?;',ic)
c_eval('SCpot = dcv?;',ic)
c_eval('zphase = zphase?;',ic)
[starttime1,endtime1,starttime3,endtime3] = mms.probe_align_times(Exyz,Bxyz,SCpot,zphase,1);
h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('irf_pl_mark(h(?),tint_waves{!},[1 0 0]);',1:numel(h),1:numel(tint_waves))

%% Calculate frequency power spectrum and wave number
probecomb = 5;
trange = tint_waves{3};
c_eval('V6 = dcv?;',ic)
c_eval('Bxyz = dmpaB?;',ic)
c_eval('zphase = zphase?;',ic)
[fkpower,freq,wavenumber] = mms.fk_powerspectrum(probecomb,trange,V6,Bxyz);

hca = subplot(1,1,1);
pcolor(hca,wavenumber,freq,fkpower);
shading(hca,'flat');
hca.XLabel.String = 'k (m^{-1})';
hca.YLabel.String = 'f (Hz)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Power/max(Power)';
hca.Title.String = [trange(1).utc ' - ' trange(2).utc];
hca.XGrid = 'on';
hca.YGrid = 'on';

if 1 % tint 3
  hold(hca,'on')
  hca.YLim = [0 120];
  xfit = [0 0.2];
  yfit = [0 40];
  vph = yfit(2)*2*pi/xfit(2);
  plot(hca,xfit,yfit,'linewidth',1)
  hold(hca,'off')
end
if 0 % tint 1
  hold(hca,'on')
  hca.YLim = [0 200];
  xfit = [0 0.2];
  yfit = [0 70];
  vph = yfit(2)*2*pi/xfit(2);
  plot(hca,xfit,yfit,'linewidth',1)
  hold(hca,'off')
end