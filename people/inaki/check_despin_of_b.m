mms.db_init('local_file_db','/Users/cecilianorgren/Data/MMS');
tint = irf.tint('2017-07-27T12:30:00.000Z/2017-07-27T14:52:00.000Z');
ic = 1:4;
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gse_srvy_l2'',tint);',ic);
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_dmpa_srvy_l2'',tint);',ic);
c_eval('zphase? = mms.db_get_variable(''mms?_ancillary_defatt'',''zphase'',tint);',ic);
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1)
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',2)
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',3)
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',4)
c_eval('gseR?_B = gseR?.resample(dmpaB?);',1:4)
[gseJcurl,divB,gseB,JxB,BdivB,gseDivPb] = c_4_j('gseR?_B','gseB?');

%%
Phase = irf.ts_scalar(EpochTT(zphase1.time),zphase1.zphase);
Data = gseJcurl.resample(Phase);
fCut = 100;
nspins = 31;
SAMPLEFREQ = Phase.time;
Out = irf_spin_epoch(Data,Phase,'NSPINS',nspins);
%Out = irf_spin_epoch(Data,Phase);
Data.name = 'In';
Out.name = 'Model';
h = irf_plot({Data,Out,Data-Out});
irf_zoom(h,'x',tint)
h(end).XTickLabelRotation = 0;
irf_legend(h(2),sprintf('nspins=%g',nspins),[0.02 0.98])