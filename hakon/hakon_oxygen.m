mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');
ic = 1:4;
tint = irf.tint('2017-08-04T10:00:03.00Z/2017-08-04T10:04:03.00Z');

c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('dobj? = dataobj(''/Volumes/Fountain/Data/MMS/mms?/fgm/brst/l2/2017/08/04/mms?_fgm_brst_l2_20170804100003_v5.98.0.cdf'');',ic)
c_eval('tt? = get_variable(dobj?,''Epoch_state'');',ic)
c_eval('rr? = get_variable(dobj?,''mms?_fgm_r_gse_brst_l2'');',ic)
c_eval('gseR? = irf.ts_vec_xyz(EpochTT(tt?.data),rr?.data(:,1:3));',ic)

%%
fs = 1/(gseB1.time(2)-gseB1.time(1));
timeline = gseB1.time(1):1/(2*fs):gseB1.time(end);
timeline = gseB1.time;
c_eval('B? = gseB?.resample(timeline);',ic)
c_eval('R? = gseR?.resample(timeline);',ic) % unnecessary

tint = gseB1.time([1 end]);
T = 10; % s
dt = T/5; % s
maxlagt = 1;
maxlags = maxlagt*fs;

t0 = gseB1.time(100);
t1 = t0;
t2 = t1 + T;
icount = 0;
while t2 < gseB1.time(end)
  icount = icount + 1
  c_eval('tmpB? = B?.tlim([t1 t2]);',ic)
  c_eval('tmpR? = R?.tlim([t1 t2]);',ic)
  
  for ic_ = 1:4
    s1 = tmpB1.x.data;
    c_eval('s2 = tmpB?.x.data;',ic_)
    %[acor,lag] = xcorr(s1-mean(s1),s2-mean(s2),maxlags,'coeff');
    [acor,lag] = finddelay_rms(s1,s2);
    [~,I] = min(abs(acor));
    timeDiff(ic_) = -lag(I)/fs;
    %plot(lag,acor,lag(I),acor(I),'*')
    %ylim([0 prctile(acor,90)])    
    %drawnow
  end  
  % Get velocity
  v_tmp = irf_4_v(tmpR1,tmpR2,tmpR3,tmpR4,t1+0.5*T+timeDiff);  
  all.tminus(icount,:) = T/2;
  all.tplus(icount,:) = T/2;
  all.v(icount,:) = v_tmp;
  all.rms(icount,:) = acor(I);
  all.dt(icount,:) = timeDiff;
  
  if 0
  h = irf_plot(2);
  hca = irf_panel('absB');
  irf_plot(hca,{tmpB1.abs,tmpB2.abs,tmpB3.abs,tmpB4.abs},'comp')
  
  hca = irf_panel('absB dt');
  irf_plot(hca,{tmpB1.abs,tmpB2.abs,tmpB3.abs,tmpB4.abs},'comp','dt',timeDiff)
  irf_zoom(h,'x',[t1 t2])
  irf_zoom(h,'y')
  drawnow
  pause(1)
  end
  t1 = t1 + dt;
  t2 = t1 + T;
  %break
  if icount>20
    %break
  end
 
end
all.t = t0+(1:icount)*dt;
gseVcorr = irf.ts_scalar(all.t,all.v);