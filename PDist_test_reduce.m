ic = 1;
tint = irf.tint('2017-07-06T16:36:00.00Z',30);
tint = irf.tint('2017-07-11T22:33:04.00Z',3);

% Load data
c_eval('iPDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0]));',ic)
c_eval('ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)
c_eval('gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)

%% Make reduced distribution
tintZoom = irf.tint('2017-07-06T16:36:05.00Z',3);
vint = [0 1000];
f1D = iPDist1.tlim(tintZoom).reduce('1D',gseB1,'vint',vint); % reduced distribution along B
lineV1 = gseVi1.dot(gseB1.norm.resample(gseVi1)); % projection of Vi on B

%% Plot
npanels = 6;
h = irf_plot(npanels);

hca = irf_panel('B');
irf_plot(hca,gseB1);
hca.YLabel.String = 'B (nT)';

hca = irf_panel('Ve');
irf_plot(hca,{gseVe1}); 
hca.YLabel.String = 'V_e (km/s)';

hca = irf_panel('Vi');
irf_plot(hca,gseVi1);
hca.YLabel.String = 'V_i (km/s)';

hca = irf_panel('iDEF');
[hout,hcb] = irf_spectrogram(hca,iPDist1.omni.deflux.specrec,'log');
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
hca.YLabel.String = {'E_i','(eV)'};   

if 1 % iPDist pa 64
  hca = irf_panel('i PA e64 deflux lowe');  
  eint = [100 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];   
end

hca = irf_panel('iLine');
irf_plot(hca,f1D.specrec('velocity_1D'));
hold(hca,'on')
irf_plot(hca,lineV1)
%irf_plot(hca,gseVi1)
hold(hca,'off')
hca.YLim = f1D.depend{1}(1,[1 end]);
hca.YLabel.String = 'v (km/s)'; 
irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  

irf_zoom(h,'x',gseB1.time)
irf_plot_axis_align
%h=irf_plot({gseB1,gseVi1,iPDist1.deflux.omni.specrec('energy'),f1D.specrec('velocity_1D')}); h(3).YScale = 'log'; %h(4).YLim = [-1000 1000];