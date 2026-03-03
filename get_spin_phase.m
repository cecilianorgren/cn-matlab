
tint = irf.tint('2017-07-27T14:16:00.00Z/2017-07-27T14:21:00.00Z');
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gse_srvy_l2'',tint);',1:4);
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4)
c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divB,gseB,JxB,gseCurvB,gseDivPb] = c_4_j('gseR?brsttime','gseB?');

spin_phase = mms.db_get_variable('mms1_ancillary_defatt','zphase',tint);
%%
tsSP = irf.ts_scalar(EpochTT((spin_phase.time)),spin_phase.zphase);

%%
figure(10); irf_plot({Jcurl,tsSP})
figure(11); plot(tsSP.resample(Jcurl).data, Jcurl.x.data,'.')
figure(12); histcn_plot([tsSP.resample(Jcurl).data, Jcurl.x.data*1e9],0:10:360,-2:0.5:20)
%%
%Phase = loadedFiles{1}.wphase(1:3000);
pha0 = 0;
phaStep = 15;
phaOffs = 0;
Phase =  [double(Phase_.time) Phase_.wphase];
aa = irf_fixed_phase_epoch(Phase,pha0,phaStep,phaOffs);