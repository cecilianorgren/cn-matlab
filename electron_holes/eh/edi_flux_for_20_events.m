% bgk_tseries_parameter_study.m
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   

ic = 1:4;
%[n_vph_all,n_phi_mult_all,n_iff_all] = size(flux_all);

doCollectCorr = 0;
if doCollectCorr
  corr_j_all = nan(n_vph_all,n_phi_mult_all,n_iff_all);
  corr_n_all = nan(n_vph_all,n_phi_mult_all,n_iff_all);
  corr_f_all = nan(n_vph_all,n_phi_mult_all,n_iff_all);
end

doPlot = 1;
if doPlot
  fig = figure(42);
  fig.Position = [100,929,1307,600];
  npanels = 4;
  h = irf_plot(npanels);
end

events = 1:20;
nEvents = numel(events);
% run through all and get correlation
for i_event = 1:nEvents
  event = events(i_event);
  sep.get_tints;
  % Load data
  c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
  c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
  c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
  c_eval('ePitch?_flux_edi = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint,?);',ic)

  c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)


  % Collect correlation data into matrix     
  % Density
      if doCollectCorr
%         % Density
%         corr_n_tmp = sum(abs(ts_n_diff.tlim(tint_phi).data).^2);
%         corr_n_all(i_vph,i_phi_mult,i_iff) = corr_n_tmp;
% 
%         % Flux
%         corr_j_tmp = sum(abs(ts_j_diff.tlim(tint_phi).data).^2);      
%         corr_j_all(i_vph,i_phi_mult,i_iff) = corr_j_tmp;
% 
%         % PSD
%         corr_f_tmp = sum(f_diff_ampl(v_ind).^2);
%         corr_f_all(i_vph,i_phi_mult,i_iff) = corr_f_tmp;            
      end
      % Averaged PSD
  vlim = [-40 -6]; % interval in which we check the difference

  % Plot TSeries
  if doPlot
    %%
    isub = 1;
    hca = h(isub); isub = isub + 1;
    irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp')
    hca.YLabel.String = 'E_{||}';
    irf_legend(hca,{sprintf('mms %g',1),sprintf('mms %g',2),sprintf('mms %g',3),sprintf('mms %g',4)},[0.02 0.98]) 

    hca = h(isub); isub = isub + 1;
    palim = 180-[11.25 0];    
    irf_plot(hca,{ePitch1_flux_edi.palim(palim),ePitch2_flux_edi.palim(palim),ePitch3_flux_edi.palim(palim),ePitch4_flux_edi.palim(palim)},'comp')
    hca.YLabel.String = {'j_{e}^{EDI}',sprintf('[%.2f,%.2f]',palim(1),palim(2))};
    irf_legend(hca,{sprintf('mms %g',1),sprintf('mms %g',2),sprintf('mms %g',3),sprintf('mms %g',4)},[0.02 0.98]) 
    irf_zoom(h,'x',tint_phi)
    
    hca = h(isub); isub = isub + 1;
    palim = [0 11.25];    
    irf_plot(hca,{ePitch1_flux_edi.palim(palim),ePitch2_flux_edi.palim(palim),ePitch3_flux_edi.palim(palim),ePitch4_flux_edi.palim(palim)},'comp')    
    hca.YLabel.String = {'j_{e}^{EDI}',sprintf('[%.2f,%.2f]',palim(1),palim(2))};
    irf_legend(hca,{sprintf('mms %g',1),sprintf('mms %g',2),sprintf('mms %g',3),sprintf('mms %g',4)},[0.02 0.98]) 
    
    
    irf_zoom(h,'x',tint_phi)

    %pause
      
  end
end
