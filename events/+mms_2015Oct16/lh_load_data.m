tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'); % magnetosphere-magnetosheath-magnetosphere

disp('Loading magnetic field...')
c_eval('dmpaB?brst=mms.db_get_ts(''mms?_dfg_brst_ql'',''mms?_dfg_brst_dmpa'',tint);');
%c_eval('dmpaB?brst.name = ''B? brst DMPA'';')

scmPath = '/Users/Cecilia/Data/MMS/2015Oct16/scm/';
scmFiles = dir([scmPath '*.mat']);
for ii = 1:numel(scmFiles)
  load([scmPath scmFiles(ii).name]);
  c_eval('dmpaB?scm = Bscm; dmpaB?scm.units = ''nT''; dmpaB?scm.name = ''scm B?''; dmpaB?scm.coordinateSystem = ''DMPA'';',ii)
end

disp('Loading electric field...')
c_eval('dslE?brst=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);');
%c_eval('dslE?3Dbrst=mms.db_get_ts(''mms?_edp_brst_ql_dce'',''mms?_edp_dce_xyz_dsl'',tint);');

%% Make E/B spectra
nfft = 512;
c_eval('fsE = 1/(dslE?brst.time(2)-dslE?brst.time(1)); pfftE?ts = irf_powerfft(dslE?brst,nfft,fsE,0.5);')
c_eval('fsB = 1/(dmpaB?scm.time(2)-dmpaB?scm.time(1)); pfftB?ts = irf_powerfft(dmpaB?scm,nfft,fsB,0.5);')

%% Make E/B spectra
nfft = 512;
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:38.00Z')+[-1 1];
c_eval('fsE = 1/(dslE?brst.time(2)-dslE?brst.time(1)); pfftE? = irf_powerfft(dslE?brst.abs.tlim(tint),nfft,fsE,0.5);')
c_eval('fsB = 1/(dmpaB?scm.time(2)-dmpaB?scm.time(1)); pfftB? = irf_powerfft(dmpaB?scm.abs.tlim(tint),nfft,fsB,0.5);')

c_eval('pfftB?.p =pfftB?.p{1};'); 
c_eval('pfftE?.p =pfftE?.p{1};'); 
c_eval('pfftE?.f_units = ''Hz''; pfftE?.f_label = ''f [Hz]''; pfftE?.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)
c_eval('pfftB?.f_units = ''Hz''; pfftB?.f_label = ''f [Hz]''; pfftB?.p_label = {''log_{10} B^2'',''(nT)^2/Hz''};',ic)

%%
nfft = 2048;
c_eval('fsE = 1/(dslE?hmfe.time(2)-dslE?hmfe.time(1)); pfftE?hmfe = irf_powerfft(dslE?hmfe.abs,nfft,fsE,0.5);')
%c_eval('fsB = 1/(dmpaB?scm.time(2)-dmpaB?scm.time(1)); pfftB?ts = irf_powerfft(dmpaB?scm,nfft,fsB,0.5);')


%%
ic = 1:4;
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:38.00Z')+[-1 1];
c_eval('wavB? = irf_wavelet(dmpaB?scm.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
c_eval('wavE? = irf_wavelet(dslE?brst.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)

c_eval('wavB?.p =wavB?.p{1};'); 
c_eval('wavE?.p =wavE?.p{1};'); 
c_eval('wavE?.f_units = ''Hz''; wavE?.f_label = ''f [Hz]''; wavE?.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)
c_eval('wavB?.f_units = ''Hz''; wavB?.f_label = ''f [Hz]''; wavB?.p_label = {''log_{10} B^2'',''(nT)^2/Hz''};',ic)
