totDir = '/data/mms/mms1/feeps/brst/l1b/electron/2015/08/28';
fileName = 'mms1_feeps_brst_l1b_electron_20150828150624_v2.2.1.cdf';
dobj=dataobj([totDir filesep fileName]);

%%
totDir = '/data/mms/mms1/feeps/brst/l1b/electron/2015/08/28';
fileName = 'mms1_feeps_brst_l1b_electron_20150828151224_v2.2.1.cdf';
dobj=dataobj([totDir filesep fileName]);
%%
totDir = '/data/mms/mms1/fpi/brst/l1b/des-moms/2015/08/15';
fileName = 'mms1_fpi_brst_l1b_des-moms_20150815122000_v0.1.1.cdf';
dobj=dataobj([totDir filesep fileName]);

%%
totDir = '/data/mms/mms1/feeps/srvy/l1b/electron/2015/08';
fileName = 'mms1_feeps_srvy_l1b_electron_20150819000000_v2.2.1.cdf';
dobj=dataobj([totDir filesep fileName]);

%%
totDir = '/data/mms/mms1/dsp/fast/l2/epsd/2015/08';
fileName = 'mms1_dsp_fast_l2_epsd_20150814_v0.6.0.cdf';
dobj=dataobj([totDir filesep fileName]);

tintEpoch = toepoch([2015 08 28 15 06 48.5; 2015 08 28 15 06 50])';
tint = EpochTT(irf_time(tintEpoch,'epoch>utc'));
epsd_omni = cn_get_ts('mms1_dsp_fast_l2_epsd','mms1_dsp_epsd_omni',tint(1));
pcolor(log10(epsd_omni.data)'); shading flat
%%
totDir = '/data/mms/mms1/fpi/brst/l1b/des-moms/2015/08/15';
fileName = 'mms1_fpi_brst_l1b_des-moms_20150815123000_v0.1.1.cdf';
dobj=dataobj([totDir filesep fileName]);
%%
totDir = '/data/mms/mms1/fpi/brst/l1b/des-dist/2015/08/15';
fileName = 'mms1_fpi_brst_l1b_des-dist_20150815123000_v0.1.1.cdf';
dobj=dataobj([totDir filesep fileName]);
%%
tint = EpochTT(irf_time(toepoch([2015 08 15 12 30 00]),'epoch>utc'));
desSkyMap = cn_get_ts('mms1_fpi_brst_l1b_des-dist','mms1_des_brstSkyMap_dist',tint(1));
%%
totDir = '/data/mms/mms1/feeps/srvy/l1b/electron/2015/08';
fileName = 'mms1_dsp_fast_l2_epsd_20150814_v0.6.0.cdf';
dobj=dataobj([totDir filesep fileName]);

%cn_get_ts;