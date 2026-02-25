%% Fujimoto 2006, how he determines the excited mode
% ion-electron instability
[wr_1,wi_1,k_1,f_1] = paper_electron_acceleration.dispersion_solver({'fujimoto_ie_quasineutral'},100,[0.08 2 100]);
[wr_2,wi_2,k_2,f_2] = paper_electron_acceleration.dispersion_solver({'fujimoto_ie_not_quasineutral'},100,[0.08 2 100]);
[n1,T1,m1,q1,vd1,vt1,wp1,Lin1,Ld1] = f_inp('fujimoto_ie_quasineutral');
[n2,T2,m2,q2,vd2,vt2,wp2,Lin2,Ld2] = f_inp('fujimoto_ie_not_quasineutral');
%%
colors = pic_colors('matlab');
h = setup_subplots(1,2);
  
  hca = h(1); % inptu distributions
  fiscale = 10;
  plot(hca,f_1{1},f_1{2},'displayname','f_{e1}','color',colors(1,:),'linestyle','-');
  hold(hca,'on')
  plot(hca,f_1{1},f_1{3}/fiscale,'displayname','f_{i1}','color',colors(1,:),'linestyle','-.');    
  plot(hca,f_2{1},f_2{2},'displayname','f_{e2}','color',colors(2,:),'linestyle','-');
  plot(hca,f_2{1},f_2{3}/fiscale,'displayname','f_{i2}','color',colors(2,:),'linestyle','-.');
  hold(hca,'off')    
  
  hca.XLabel.String = 'v_ (km/s)';
  hca.YLabel.String = 'f (s/m^4)';
  hca.Title.String = 'input distributions';
  hleg = legend(hca);
  hleg.Title.String = {'Species,','e - electrons','i - ions'};
  hleg.Box = 'off'; 
  
  hca = h(2);
  plot(hca,k_1*1e3,wr_1/(2*pi)*1e-3,'k','displayname','f_r','color',colors(1,:))
  hold(hca,'on')
  plot(hca,k_1*1e3,wi_1/(2*pi)*1e-3,'k--','displayname','f_i','color',colors(1,:))
  plot(hca,k_2*1e3,wr_2/(2*pi)*1e-3,'k','displayname','f_r','color',colors(2,:))
  plot(hca,k_2*1e3,wi_2/(2*pi)*1e-3,'k--','displayname','f_i','color',colors(2,:))
  hold(hca,'off')
 
  hca.XLabel.String = 'k (km^{-1})';
  hca.YLabel.String = 'f (kHz)';  
  hca.Title.String = 'dispersion relation';
  hca.YLim = [0 0.2];
  hca.XLim = [0 4];
  hleg = legend(hca);
  hleg.Box = 'off';
  hleg.Title.String = {'Complex frequencies','r - real','i - imganiary'};
 