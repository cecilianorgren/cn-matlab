%% General
% Access overview plots made by Wenya Li
% on spis: /data/mms/irfu/plots/summary_plot/

% Mount remote disc locally to access database
% I do:
% sshfs cecilia@spis.irfu.se:/data/mms/ /Users/cecilia/Discs/spis -o uid=501,Compression=no,volname=spis,negative_vncache
% where /Users/cecilia/Discs/spis is a local empty folder I set up and 
% /data/mms/ is the local file path to the mms data on spis
% Locally it exists as /Users/cecilia/Discs/spis_data_mms

% It would probably be faster to run Matlab on spis

%% Initialize database
pathToDatabse = '/Users/cecilia/Data/MMS'; % locally stored files
pathToDatabse = '/Users/cecilia/Discs/spis'; % files stored on spis, which I mounted locally using sshfs
mms.db_init('local_file_db',pathToDatabse);
db_info = datastore('mms_db');

%% Specify time interval in which to search for existing files
tint = irf.tint('2017-06-01T00:00:00.00Z/2017-07-01T00:00:00.00Z');

%% Make a list of all existing files in the database
% Use some part of filename that should exist for all times, for example
% magnetic field.
% Example of full name: mms1_fgm_brst_l2_20170706135303_v5.92.0.cdf
% 20170706135303 - yyyymmddHHMMSS - specifies the burst interval

% Search for EDI file name instead, because we're only intersted in those
% intervals. There can be two EDI name types I think.
tic % Time how log time it takes to make file list
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);
toc

%% The above takes forever, not sure why...
% Attempt at faster version.
disp('Compiling file list.')
tic;
%dd = dir([db_info.local_file_db_root filesep 'mms1' filesep 'fgm' filesep 'brst' filesep '**' filesep '*mms1_fgm_brst_l2*.cdf']);
dd = dir([db_info.local_file_db_root filesep 'mms1' filesep 'edi' filesep 'brst' filesep '**' filesep '*mms1_edi_brst_l2_amb*.cdf']);
toc
%% This may contain two file versions of the same time interval 
% (e.g. v5.1.0 and v5.1.1), but we only want the latest version. So we
% remove duplicate time stamps.
fileNamesCell = {dd.name}';
%fileIdsCell = cellfun(@(x) x(18:31),fileNamesCell,'UniformOutput',false);
%fileIdsString = unique(string(fileIdsCell));
tmp = cellfun(@(x) regexp(x, '_(\d*)_v','tokens'),fileNamesCell,'UniformOutput',false);
fileIdsCell = cellfun(@(x) x{1}, tmp,'UniformOutput',false);
fileIdsString = unique(string(fileIdsCell));

%% Load data
fileId = 1; % can be adapted into loop

tint_str = sprintf('%s-%s-%sT%s:%s:%s',fileIdsString{fileId}(1:4),fileIdsString{fileId}(5:6),fileIdsString{fileId}(7:8),fileIdsString{fileId}(9:10),fileIdsString{fileId}(11:12),fileIdsString{fileId}(13:14));
files = mms.db_list_files('mms1_fgm_brst_l2',EpochTT(tint_str) + [0 1]);
tint = [files(1).start files(1).stop]; % the time interval of that file id
%%
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');

ic = 1; % specify what spacecraft you want to load data from

% c_eval replaces every ? with the number in the 2nd argument, and any !
% with an optional third argument, and then evaluates the statement. Any 
% single ' must be replaced by two ''. Test 
% >> c_eval('disp(''!?'')',1:2,4:5)

disp('Loading data.')
% Magnetic field
%c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
%c_eval('gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);

% Electric field
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
%c_eval('dslE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
%c_eval('gsmE? = c_coord_trans(''GSE'',''GSM'',gseE?);',ic)

% Spacecraft potential
% The spacecraft potential can be used as a proxy for the plasma density
% (https://www.cluster.irfu.se/efw/ops/dummies/index.html) and has a higher
% sampling rate than the 'plasma instrument density'.
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

% Density
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);

% Velocity
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);
%c_eval('gsmVe? = c_coord_trans(''GSE'',''GSM'',gseVe?);',ic)
%c_eval('gsmVi? = c_coord_trans(''GSE'',''GSM'',gseVi?);',ic)

% EDI
c_eval('ePitch?_flux_edi            = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint,?);',ic)
c_eval('ePitch?_flux_edi_err        = mms.get_data(''Flux-err-amb-pm2_edi_brst_l2'',tint,?);',ic)
c_eval('ePitch?_flux_edi_err_plus   = ePitch?_flux_edi + (ePitch?_flux_edi_err*+1);',ic)
c_eval('ePitch?_flux_edi_err_minus  = ePitch?_flux_edi + (ePitch?_flux_edi_err*-1);',ic)        

%c_eval('ePitch?_counts_edi = mms.get_data(''Counts-amb-pm2_edi_brst_l1a'',tint,?);',ic)

% Spacecraft position
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4)

disp('Done loading data.')

%% Derive quanitites 
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

% Integrate electric field
c_eval('intE?par = irf_integrate(gseE?par);',ic)

% high-pass filter integrated electric field
flow = 100; % Hz
fhigh = 0; % Hz, 0 = Inf
fsample = []; % if empty, calculated automatically
ffiltorder = 5;
c_eval('filtintE?par = irf_filt(intE?par,flow,fhigh,fsample,ffiltorder);',ic) % help irf_filt

% Extract EDI flux at 0 and 180 pitch angles
c_eval('jedi_apar = ePitch?_flux_edi.palim([168.75 180]);',ic)
c_eval('jedi_par = ePitch?_flux_edi.palim([0 11.25]);',ic)

%% Make wave power spectrogram, for easy overview of wiggles  
nfft = 256; % how many data points one power point should include
sfreqE = 1/(gseE1.time(2)-gseE1.time(1)); % Hz
sfreqFlux = 1/(ePitch1_flux_edi.time(2)-ePitch1_flux_edi.time(1)); % Hz
overlap = 0.50;
c_eval('fftE? = irf_powerfft(intE?par,nfft,sfreqE,[overlap]);?',ic)
c_eval('fftFlux?_0 = irf_powerfft(ePitch?_flux_edi.palim(0),nfft,sfreqFlux,[overlap]);?',ic)
c_eval('fftFlux?_180 = irf_powerfft(ePitch?_flux_edi.palim(180),nfft,sfreqFlux,[overlap]);?',ic)
  
%% Cross correlation of data
% try matlab function xcorr to correlate filtintE?par and jedi_par/jedi_apar
% to correlate with xcorr, the two data need the same number of data points. 
% need to downsample potential to jEDI timeline.
% E.g.
% A = filtintE1par.resample(ePitch1_flux_edi);
% C = xcorr(A,B,0);

% Something with wave spectra cross-correlation...
tint_zoom = tint; irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.620Z');
A = jedi_apar.tlim(tint_zoom).data;
B = filtintE1par.resample(jedi_apar).tlim(tint_zoom).data;
%pxy = cpsd(A,B);
%tsC = irf.ts_scalar(jedi_apar.time,[real(pxy), imag(pxy)]);
[S,F,T] = xspectrogram(A,B,100);
hca = subplot(4,1,1);
pcolor(hca,T,F,log10(S))
shading flat
hca = subplot(4,1,2);
plot(hca,A)
hca = subplot(4,1,3);
plot(hca,B)

hca = subplot(4,1,4);
plot(hca,gseE1par.resample(jedi_apar).tlim(tint).data)

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
linkprop(h(:),{'XLim'});
%% Other useful function
% TSeries.tlim(tint)

%% Make figure
% What spacecraft to plot data for
ic = 1;
% Initiate figure with a number of panels
npanels = 10;
h = irf_plot(npanels); 


if 1 % B GSE
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % scPot - spacecraft potential, proxy for -density
  hca = irf_panel('scPot');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{-1*scPot?},''comp'');',ic)
  hca.YLabel.String = {'V_{SC}','(V)'};
end
if 1 % Vi  
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % E par
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 1 % Electron flux from EDI instrument, parallel and antiparallel to B
  hca = irf_panel('j EDI par apar');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{jedi_par,jedi_apar},''comp'');',ic)
  hca.YLabel.String = {'j_{EDI}','(cm^{-2} s^{-1})'};
end
if 1 % E power spectrogram
  hca = irf_panel('Epar fft');
  c_eval('irf_spectrogram(hca,fftE?.t,fftE?.p{1},fftE?.f)',ic);
  hca.YScale = 'log';
  hc = colorbar('peer',hca);
  hc.Label.String = '(mV/m)^2/Hz';    
  %hca.YTick = [1 10 100 1000 10000];
  hca.YLabel.String = 'f [Hz]';
end
if 1 % j EDI par power spectrogram
  hca = irf_panel('j EDI par fft');
  c_eval('irf_spectrogram(hca,fftFlux?_0.t,fftFlux?_0.p{1},fftFlux?_0.f)',ic);
  hca.YScale = 'log';
  hc = colorbar('peer',hca);
  hc.Label.String = '(...)^2/Hz';    
  %hca.YTick = [1 10 100 1000 10000];
  hca.YLabel.String = 'f [Hz]';
end
if 1 % j EDI apar power spectrogram
  hca = irf_panel('j EDI apar fft');
  c_eval('irf_spectrogram(hca,fftFlux?_180.t,fftFlux?_180.p{1},fftFlux?_180.f)',ic);
  hca.YScale = 'log';
  hc = colorbar('peer',hca);
  hc.Label.String = '(...)^2/Hz';    
  %hca.YTick = [1 10 100 1000 10000];
  hca.YLabel.String = 'f [Hz]';
end

irf_colormap('waterfall')
irf_plot_axis_align
irf_zoom(h,'x',tint)