%% Intialize database
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
db_info = datastore('mms_db');   

%% Time interval 1
tint = irf.tint('2017-08-20T01:52:43.00Z/2017-08-20T01:53:43.00Z');
eis_omni = mms.get_data('Omnifluxproton_epd_eis_brst_l2',tint,3);

h = irf_plot(1);
hca = irf_panel('eis omni');
[hout,hcb] = irf_spectrogram(hca,eis_omni.specrec,'log');

%% Time interval 2
tint = irf.tint('2017-06-11T17:24:53.00Z/2017-06-11T17:25:53.00Z');
eis_omni = mms.get_data('Omnifluxproton_epd_eis_brst_l2',tint,3);

h = irf_plot(1);
hca = irf_panel('eis omni');
[hout,hcb] = irf_spectrogram(hca,eis_omni.specrec,'log');

%% Time interval 3
tint = irf.tint('2020-06-26T00:00:00.00Z/2020-06-26T02:00:00.00Z');
eis_omni = mms.get_data('Omnifluxproton_epd_eis_srvy_l2',tint,1);

h = irf_plot(1);
hca = irf_panel('eis omni');
[hout,hcb] = irf_spectrogram(hca,eis_omni.specrec,'log');

%% Time interval 3, brst
tint = irf.tint('2020-06-26T00:00:00.00Z/2020-06-26T02:00:00.00Z');
eis_omni = mms.get_data('Omnifluxproton_epd_eis_brst_l2',tint,1);

h = irf_plot(1);
hca = irf_panel('eis omni');
[hout,hcb] = irf_spectrogram(hca,eis_omni.specrec,'log');