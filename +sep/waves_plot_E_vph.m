%% Read existing data
file_list = dir([pathData '/*.txt']);
table = [];
n_events = numel(file_list);

time_center = '';
time_start = '';
time_stop = '';
vph = nan(n_events,3);
phi = nan(n_events,1);
length = nan(n_events,1);
correlation  = nan(n_events,1);

for i_event = 1:n_events
  fid = fopen([pathData file_list(i_event).name],'r');  
  tmp_data = textscan(fid,data_format_read);
  %time_start_tmp = tmp_data{3}{1};
  %time_stop_tmp = tmp_data{4}{1};
  %time_center_tmp = tmp_data{5}{1};
  time_start(i_event,:) = tmp_data{3}{1};
  time_stop(i_event,:) = tmp_data{4}{1};
  time_center(i_event,:) = tmp_data{5}{1};
  vph(i_event,1:3) = [tmp_data{6} tmp_data{7} tmp_data{8}];
  phi(i_event,1:4) = [tmp_data{10} tmp_data{11} tmp_data{12} tmp_data{13}];
  length(i_event,1:4) = [tmp_data{14} tmp_data{15} tmp_data{16} tmp_data{17}];
  correlation(i_event,1:4) = [tmp_data{18} tmp_data{19} tmp_data{20} tmp_data{21}];
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
%% Make TSeries of extisting data
tsVph = irf.ts_vec_xyz(time_center,vph);
tsPhi = irf.ts_scalar(time_center,phi);
tsLength = irf.ts_scalar(time_center,length);
tsCorrelation = irf.ts_scalar(time_center,correlation);

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

hmark_phi = irf_pl_mark(h(1),[tint_phi(1).epochUnix tint_phi(2).epochUnix]);


if ~isempty(time_start)
  hmark = irf_pl_mark(h(1),[time_start.epochUnix time_stop.epochUnix],'b');
end