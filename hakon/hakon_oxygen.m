mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');
ic = 1:4;
tint = irf.tint('2017-08-04T10:00:03.00Z/2017-08-04T10:04:03.00Z');
tint = irf.tint('2018-07-05T20:19:00.00Z/2018-07-05T20:23:45.00Z');
tint = irf.tint('2018-07-05T00:00:00.00Z/2018-07-05T24:00:00.00Z');
tint = irf.tint('2018-07-05T20:23:00.00Z/2018-07-05T20:23:45.00Z');

tint = irf.tint('2018-07-05T20:21:00.00Z/2018-07-05T20:24:00.00Z');

c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
c_eval('nO? = mms.get_data(''Noplus_hpca_brst_l2'',tint,?);',ic);

%c_eval('dobj? = dataobj(''/Volumes/Fountain/Data/MMS/mms?/fgm/brst/l2/2017/08/04/mms?_fgm_brst_l2_20170804100003_v5.98.0.cdf'');',ic)
%c_eval('tt? = get_variable(dobj?,''Epoch_state'');',ic)
%c_eval('rr? = get_variable(dobj?,''mms?_fgm_r_gse_brst_l2'');',ic)
%c_eval('gseR? = irf.ts_vec_xyz(EpochTT(tt?.data),rr?.data(:,1:3));',ic)
R = mms.get_data('R_gse',tint);
if size(R.gseR1,2) == 4
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre  
else
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
end

c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s

c_eval('tic; [iPDist?,iPDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)
%c_eval('iPDist?.data(iPDist?.data<iPDistErr?.data*1.1) = 0;',ic)
tic; if1Dz = iPDist1.reduce('1D',[1 0 100]); toc % reduced distribution along B

%% Check data
dobj = dataobj('/Volumes/Fountain/Data/MMS/mms1/fgm/brst/l2/2018/07/05/mms*.cdf');
tt = get_variable(dobj,'Epoch');
gsmb = get_variable(dobj,'mms1_fgm_b_gsm_brst_l2');
EpochTT(tt.data)
gsmB1 = irf.ts_vec_xyz(EpochTT(tt.data),gsmb.data(:,1:3));

%%
tint = irf.tint('2018-07-05T20:23:00.00Z/2018-07-05T20:23:45.00Z');
[out,eigenVal,eigenVec] = irf_minvar(gsmB1.tlim(tint),'td'); % Same as the interactive w/ column 3
lmn = [eigenVec(1,:);eigenVec(2,:);eigenVec(3,:)];
gsmB1lmn = gsmB1*lmn';

%% RMS correlation
fs = 1/(gseB1.time(2)-gseB1.time(1));
timeline = gseB1.time(1):1/(2*fs):gseB1.time(end);
timeline = gseB1.time;
c_eval('B? = gseB?.resample(timeline);',ic)
c_eval('R? = gseR?.resample(timeline);',ic) % unnecessary

tint = gseB1.time([1 end]);
T = 5; % s
dt = T/2; % s, window overlap
maxlagt = 1;
maxlags = fix(maxlagt*fs);

clear all_data

t0 = gseB1.time(100);
t1 = t0;
t2 = t1 + T;
icount = 0;
%fprintf('\n it/nt = %g,',icount)
while t2 < gseB1.time(end)
  icount = icount + 1;
  fprintf('%g,',icount)
  c_eval('tmpB? = B?.tlim([t1 t2]);',ic)
  c_eval('tmpR? = R?.tlim([t1 t2]);',ic)
  
  for ic_ = 1:4
    s1 = tmpB1.abs.data;
    c_eval('s2 = tmpB?.abs.data;',ic_)
    %[acor,lag] = xcorr(s1-mean(s1),s2-mean(s2),maxlags,'coeff');
    [acor,lag] = finddelay_rms(s1,s2,maxlags);
    [~,I] = min(abs(acor));
    timeDiff(ic_) = -lag(I)/fs;
    %plot(lag,acor,lag(I),acor(I),'*')
    %ylim([0 prctile(acor,90)])    
    %drawnow
  end  
  % Get velocity
  v_tmp = irf_4_v(tmpR1,tmpR2,tmpR3,tmpR4,t1+0.5*T+timeDiff);  
  v_norm = v_tmp/norm(v_tmp);
  all_data.tminus(icount,:) = T/2;
  all_data.tplus(icount,:) = T/2;
  all_data.v(icount,:) = v_tmp;
  all_data.rms(icount,:) = acor(I);
  all_data.dt(icount,:) = timeDiff;
  
  if 1 % Plot results
    h = irf_plot(2);

    hca = irf_panel('absB');
    irf_plot(hca,{tmpB1.abs,tmpB2.abs,tmpB3.abs,tmpB4.abs},'comp')

    hca = irf_panel('absB dt');
    irf_plot(hca,{tmpB1.abs,tmpB2.abs,tmpB3.abs,tmpB4.abs},'comp','dt',timeDiff)
    leg{1,1} = sprintf('dt = [%.0f,%.0f,%.0f,%.0f] ms',timeDiff(1)*1e3,timeDiff(2)*1e3,timeDiff(3)*1e3,timeDiff(4)*1e3);
    leg{2,1} = sprintf('v = %.0f x [%.2f,%.2f,%.2f] km/s',norm(v_tmp),v_norm(1),v_norm(2),v_norm(3));
    irf_legend(hca,leg,[0.98 0.98])

    irf_zoom(h,'x',[t1 t2])
    irf_zoom(h,'y')
    drawnow
    pause
  end
  t1 = t1 + dt;
  t2 = t1 + T;
  %break
  if icount>20
    %break
  end
 
end
disp('Done.')
all_data.t = t0+(1:icount)*dt+0.5*T;
gseV = irf.ts_scalar(all_data.t,all_data.v); gseV.name = 'V_{corr}';
gseDt = irf.ts_scalar(all_data.t,all_data.dt); gseDt.name = 'dt_{corr}';
gseCorr = irf.ts_scalar(all_data.t,all_data.rms); gseCorr.name = 'C_{corr}';

%% PLot results of correlation
npanels = 8;
h = irf_plot(npanels);

if 1 % B
  hca = irf_panel('gseB1');
  irf_plot(hca,{gseB1});
  hca.YLabel.String = 'B (nT)';
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 1 % E
  hca = irf_panel('gseE1');
  irf_plot(hca,{gseE1});
  hca.YLabel.String = 'E (mV/m)';
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 1 % n
  hca = irf_panel('n');
  nOmult = 10;
  irf_plot(hca,{ni1,ne1,nO1*nOmult},'comp');
  hca.YLabel.String = 'n (cm^{-3})';
  irf_legend(hca,{'n_{i,fpi}','n_{e,fip}',sprintf('%gn_{O+,hpca}',nOmult)},[0.02 0.98])
end
if 1 % iDEF omni
  hca = irf_panel('iDEF');  
  [hout,hcb] = irf_spectrogram(hca,iPDist1.deflux.omni.specrec,'log');
%   hold(hca,'on')
%   lineScpot = irf_plot(hca,scPot1,'k');
%   lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
%   hold(hca,'off')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_i','(eV)'};   
  colormap(hca,'jet')
end
if 1 % i psd vz
  hca = irf_panel('iFred z');
  irf_spectrogram(hca,if1Dz.specrec('velocity_1D'));
  hold(hca,'on')
  irf_plot(hca,gseVExB1.resample(if1Dz).z,'linewidth',0.5)
  %irf_plot(hca,gseVExB1.resample(if1Dz).z*sqrt(16),'linewidth',0.5)
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hca.YLim = if1Dz.depend{1}(1,[1 end]);
  hca.YLabel.String = 'v_{iz} (km/s)'; 
  colormap(hca,'jet')
  irf_legend(hca,{'v_{ExB}'},[0.98 0.98])
end
if 1 % dt
  hca = irf_panel('dt');
  irf_plot(hca,{gseDt})
  irf_legend(hca,{'1-1','2-1','3-1','4-1'},[0.98 0.98])
end
if 1 % Corr
  hca = irf_panel('Corr');
  irf_plot(hca,{gseCorr})
  %irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 1 % V
  hca = irf_panel('V');
  irf_plot(hca,{gseV})
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end

irf_zoom(h,'x',tint)
irf_plot_axis_align

for ip = 1:npanels
  h(ip).YLabel.Interpreter = 'tex';
end