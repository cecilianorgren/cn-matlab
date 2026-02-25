%% Read existing data
file_list = dir([pathData '/*.txt']);
table = [];
n_events = numel(file_list);

time_start_vic = '';
time_stop_vic = '';
time_center = '';
time_start = '';
time_stop = '';
time_ref = '';
vph_all = nan(n_events,6);
vph_av = nan(n_events,1);
vph_4xcorr = nan(n_events,1);
vph_4xcorr_par = nan(n_events,3);
phi = nan(n_events,1);
length = nan(n_events,1);
correlation  = nan(n_events,1);
avB = nan(n_events,3);
sc_obs = nan(n_events,4);


for i_event = 1:n_events
  fid = fopen([pathData file_list(i_event).name],'r');  
  tmp_data_orig = textscan(fid,data_format_read);
  tmp_data = tmp_data_orig;
  
  time_start_vic(i_event,:) = tmp_data{1}{1}; tmp_data = tmp_data(2:end);
  time_stop_vic(i_event,:) = tmp_data{1}{1}; tmp_data = tmp_data(2:end);
  time_start(i_event,:) = tmp_data{1}{1}; tmp_data = tmp_data(2:end);
  time_stop(i_event,:) = tmp_data{1}{1}; tmp_data = tmp_data(2:end);
  time_center(i_event,:) = tmp_data{1}{1}; tmp_data = tmp_data(2:end);
  tmp_data = tmp_data(5:end); % tref    
  avB(i_event,1:3) = [tmp_data{1} tmp_data{2} tmp_data{3}]; tmp_data = tmp_data(4:end);
  sc_obs(i_event,1:4) = [tmp_data{1} tmp_data{2} tmp_data{3} tmp_data{4}]; tmp_data = tmp_data(5:end);  
  C_all(i_event,1:6) = [tmp_data{1} tmp_data{2} tmp_data{3} tmp_data{4} tmp_data{5} tmp_data{6}]; tmp_data = tmp_data(7:end);  
  dt_all(i_event,1:6) = [tmp_data{1} tmp_data{2} tmp_data{3} tmp_data{4} tmp_data{5} tmp_data{6}]; tmp_data = tmp_data(7:end);  
  dRpar_all(i_event,1:6) = [tmp_data{1} tmp_data{2} tmp_data{3} tmp_data{4} tmp_data{5} tmp_data{6}]; tmp_data = tmp_data(7:end);  
  vph_all(i_event,1:6) = [tmp_data{1} tmp_data{2} tmp_data{3} tmp_data{4} tmp_data{5} tmp_data{6}]; tmp_data = tmp_data(7:end);  
  vph_av(i_event,1) = [tmp_data{1}]; tmp_data = tmp_data(2:end);  
  vph_4xcorr(i_event,1:3) = [tmp_data{1} tmp_data{2} tmp_data{3}]; tmp_data = tmp_data(4:end);  
  vph_4xcorr_par(i_event,1) = [tmp_data{1}]; tmp_data = tmp_data(2:end);  
  phi(i_event,1:4) = [tmp_data{1} tmp_data{2} tmp_data{3} tmp_data{4}]; tmp_data = tmp_data(5:end);    
  length(i_event,1:4) = [tmp_data{1} tmp_data{2} tmp_data{3} tmp_data{4}]; tmp_data = tmp_data(5:end);  
  
  fclose(fid);
end

if ~isempty(time_start)
  time_start = EpochTT(time_start);
  time_stop = EpochTT(time_stop);
  time_center = EpochTT(time_center);
else
  time_start = EpochTT([]);
  time_stop = EpochTT([]);
  time_center = EpochTT([]);
end

phi(sc_obs==0) = NaN;
%% Make TSeries of extisting data
tsVph = irf.ts_scalar(time_center,vph_av);
tsPhi = irf.ts_scalar(time_center,phi);
tsLength = irf.ts_scalar(time_center,length);
%tsCorrelation = irf.ts_scalar(time_center,correlation);

%%
figure(75);
npanels = 5;
h = irf_plot(npanels);

hca = irf_panel('epar');
irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp');
hca.YLabel.String = 'E_{||}';

hca = irf_panel('eperp');
irf_plot(hca,{gseE1perp.abs,gseE2perp.abs,gseE3perp.abs,gseE4perp.abs},'comp');
hca.YLabel.String = '|E_{\perp}|';

hca = irf_panel('vph');
irf_plot(hca,{tsVph});
hca.YLabel.String = 'v_{ph}';

hca = irf_panel('phi');
irf_plot(hca,{tsLength});
hca.YLabel.String = 'phi';

hca = irf_panel('corr');
irf_plot(hca,{tsCorrelation});
hca.YLabel.String = 'Correlation';

hmark_phi = irf_pl_mark(h(1),[tint_phi(1).epochUnix tint_phi(2).epochUnix]); hmark_phi.FaceAlpha = 0.5;


if ~isempty(time_start)
  hmark = irf_pl_mark(h(1),[time_start.epochUnix time_stop.epochUnix],'b');
end
irf_zoom(h,'x',gseE1.time([1 end]))
