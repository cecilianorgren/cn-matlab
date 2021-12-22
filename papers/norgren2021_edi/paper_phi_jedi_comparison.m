% paper_figure, phi,edi comparison only
%% Load data
%ic = 1:4;
units = irf_units;
tint_brst = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tint = tint_brst + [+5 -5]; % using the above edges causes problem with new EDI files because they have different versions that adjoining file
tint_figure = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.630Z');
tint_figure = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.620Z');
tint = tint_figure + [-5 5];
localuser = datastore('local','user');
db_info = datastore('mms_db');   

c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); ',1:4);

c_eval('ePitch?_flux_edi = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint,?);',1:4)

%% Prepare data
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',1:4)

n0 = 0.04; % for normalization of density perturbation later

vph = -8500e3;
c_eval('[phi?,phi_progressive?,phi_ancillary?] = get_phi(gseE?par,vph,tint_zoom,tint_zoom);',1:4)


%% Plot, only 1 panel of EDI phi comparison
h = irf_plot(1); 
colors = mms_colors('12');
colors = [0 0 0;  0 0.4470 0.7410];
%colors = [0 0 0;  0.4660    0.6740    0.1880];
%colors = [0 0.4470 0.7410;  0.4660    0.6740    0.1880];
ylim1 = [0 3.6];
ylim2 = [0 360];
if 1 % edi, phi comparison for one spacecraft
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi phi comp');  
  set(hca,'ColorOrder',colors)
  c_eval('irf_plot(hca,{ePitch?_flux_edi.palim(palim)*1e-6},''color'',colors(1,:));',ic)
  %irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6});
  %hca.YLabel.String = {'j_e^{EDI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  %hca.YLabel.String = {'j_e^{EDI}','\theta = [168.75 180]^o','(10^6 s^{-1}cm^{-2}sr^{-1})'};  
  hca.YLabel.String = {'j_e^{EDI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};  
  hca.YLabel.Interpreter = 'tex';
  
  set(hca,'ColorOrder',colors)
 
  hca.YLim = ylim1; 
  yyaxis('right')
  ax = gca;
  irf_plot(ax,phi1)
  ax.YLabel.String = {'\phi (V)'};
  ax.YLabel.Interpreter = 'tex';
  ax.YLim = ylim2;
  %ax.YLabel.Position(1) = 1.07;
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.02 0.99],'fontsize',12,'color',[0 0 0]);
  %hca.Title.String = 'mms 1';
end
grid off
hca.Position(2) = 0.22;
hca.Position(4) = 0.65;
irf_zoom(h,'x',phi1.time)
irf_zoom(h,'x',tint_figure)