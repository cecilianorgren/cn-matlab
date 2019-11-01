% bgk_tseries_parameter_study.m
% run bgk_tseries to get some basic stuff (e.g. ts_edi_flux180/f_scale), and/or then
% load('/Users/cno062/MATLAB/cn-matlab/liouville/liouville_eh_abel_parameter_study.mat')

%flux_all{i_vph,i_phi_mult,i_iff} = fluxModel180/f_scale;
%psd_all{i_vph,i_phi_mult,i_iff} = {{v_vec*v_scale*1e-3,mod_f0},{v_vec*v_scale*1e-3,mod_f_average},{v_fpi*v_scale,1*f_fpi*1e0}};
%n_all{i_vph,i_phi_mult,i_iff} = {-1*tsDnFromPhi/nscale,-1*(tsDnModel-ntot*1e-6)/nscale};



[n_vph_all,n_phi_mult_all,n_iff_all] = size(flux_all);

doCollectCorr = 0;
if doCollectCorr
  corr_j_all = nan(n_vph_all,n_phi_mult_all,n_iff_all);
  corr_n_all = nan(n_vph_all,n_phi_mult_all,n_iff_all);
  corr_f_all = nan(n_vph_all,n_phi_mult_all,n_iff_all);
end

doPlot = 1;
if doPlot
  fig = figure(42);
  fig.Position = [1000,929,1307,600];
  npanels = 4;
  [h,h2] = initialize_combined_plot(npanels,2,1,0.6,'vertical'); % horizontal
end
      
% run through all and get correlation
for i_vph = 1:n_vph_all
  for i_phi_mult = 1:n_phi_mult_all
    for i_iff = 8%1:n_iff_all %3%[3 8] % 
      % Pick out data from cell matrix, just for ease of reading code.
      ts_n_mod = n_all{i_vph,i_phi_mult,i_iff}{2};
      ts_n_obs = n_all{i_vph,i_phi_mult,i_iff}{1};
      ts_j_mod = flux_all{i_vph,i_phi_mult,i_iff}.resample(ts_edi_flux180);
      ts_j_obs = ts_edi_flux180/f_scale;
      
      common_v_vec = psd_all{i_vph,i_phi_mult,i_iff}{3}{1}; %  from observations, fewer points
      f0_mod_resamp =  interp1(psd_all{i_vph,i_phi_mult,i_iff}{1}{1},psd_all{i_vph,i_phi_mult,i_iff}{1}{2},common_v_vec);
      fav_mod_resamp = interp1(psd_all{i_vph,i_phi_mult,i_iff}{2}{1},psd_all{i_vph,i_phi_mult,i_iff}{2}{2},common_v_vec);
      f_diff_ampl = abs(fav_mod_resamp-mean(psd_all{i_vph,i_phi_mult,i_iff}{3}{2},1)).^1;      
      v_incl = [-40 -6]; % range of velocities to include in comparison
      v_ind = intersect(find(common_v_vec>=v_incl(1)),find(common_v_vec<=v_incl(2))); % including boundaries
      
      
      
      % Calculate difference between modela and observations, TSeries and
      % correlation amplitudes.
      ts_n_diff = n_all{i_vph,i_phi_mult,i_iff}{1}-n_all{i_vph,i_phi_mult,i_iff}{2};
      ts_j_diff = flux_all{i_vph,i_phi_mult,i_iff}.resample(ts_edi_flux180); % two steps due to one being TSeries and on PDist
      ts_j_diff.data = ts_j_diff.data-ts_edi_flux180.data/f_scale;
      ts_j_diff = ts_j_diff.tlim(tint_phi);
      
      % Collect correlation data into matrix     
      % Density
      if doCollectCorr
        corr_n_tmp = sum(abs(ts_n_diff.tlim(tint_phi).data).^2);
        corr_n_all(i_vph,i_phi_mult,i_iff) = corr_n_tmp;

        % Flux
        corr_j_tmp = sum(abs(ts_j_diff.tlim(tint_phi).data).^2);      
        corr_j_all(i_vph,i_phi_mult,i_iff) = corr_j_tmp;

        % PSD
        corr_f_tmp = sum(f_diff_ampl(v_ind).^2);
        corr_f_all(i_vph,i_phi_mult,i_iff) = corr_f_tmp;            
      end
      % Averaged PSD
      vlim = [-40 -6]; % interval in which we check the difference
      
      % Plot TSeries
      if doPlot
        isub = 1;
        hca = h(isub); isub = isub + 1; hold(hca,'off')      
        irf_plot(hca,{ts_n_obs,ts_n_mod},'comp')
        hca.YLim = [-1 2];

        hca = h(isub); isub = isub + 1; hold(hca,'off')      
        irf_plot(hca,{abs(ts_n_diff).^2},'comp')
        hca.YLim = [0 2];
        irf_legend(hca,{['\Sigma_t' sprintf('|n^{mod}-n^{edi}|^2 = %g',corr_n_tmp)]},[0.02 0.98])

        hca = h(isub); isub = isub + 1; hold(hca,'off')
        irf_plot(hca,{ts_j_obs,ts_j_mod},'comp')
        hca.YLim = [0 3.5];

        hca = h(isub); isub = isub + 1; hold(hca,'off')      
        irf_plot(hca,{abs(ts_j_diff).^2},'comp')
        hca.YLim = [0 4];
        irf_legend(hca,{['\Sigma_t' sprintf('|j^{mod}-j^{edi}|^2 = %g',corr_j_tmp)]},[0.02 0.98])

        irf_zoom(h,'x',tint_phi)
        
        h(1).Title.String = sprintf('iff = %.0f, v_{ph} = %.0f km/s,  phi multiplier = %.1f',iff_all(i_iff),vph_all(i_vph)*1e-3,phi_mult_all(i_phi_mult));

        % Plot average distributions
        isub = 1;
        hca = h2(isub); isub = isub + 1; hold(hca,'off')      
        plot(hca,common_v_vec,f0_mod_resamp,...
                 common_v_vec,fav_mod_resamp,'linewidth',1.5)
        hold(hca,'on')
        plot(hca,psd_all{i_vph,i_phi_mult,i_iff}{3}{1},psd_all{i_vph,i_phi_mult,i_iff}{3}{2},'--')
        plot(hca,psd_all{i_vph,i_phi_mult,i_iff}{3}{1},mean(psd_all{i_vph,i_phi_mult,i_iff}{3}{2},1),'linewidth',1.5)
        hold(hca,'off')
        hca.XLabel.String = 'v_{||} (10^3 km/s)';
        hca.YLabel.String = 'f (s/m^4)';      
        hca.XLim = 40*[-1 1];
        hca.YLim = [0 3.5*1e-3];
        legend(hca,{'f_0','<f^{mod}>','f^{fpi}','f^{fpi}','f^{fpi}','f^{fpi}','<f^{fpi}>'})

        % Difference
        hca = h2(isub); isub = isub + 1; hold(hca,'off')      
        plot(hca,common_v_vec,f_diff_ampl,'linewidth',1.5)
        hold(hca,'on')
        
        patch_x = [common_v_vec(v_ind) common_v_vec(v_ind(end:-1:1))];
        patch_y = [f_diff_ampl(v_ind) 0*f_diff_ampl(v_ind(end:-1:1))];
        patch(patch_x,patch_y,patch_x*0)
        hold(hca,'off')
        hca.XLabel.String = 'v_{||} (10^3 km/s)';
        hca.YLabel.String = '<f^{mod}>-<f^{fpi}> (s/m^4)';        
        hca.XLim = 40*[ -1 1];
        hca.YLim = [0 3.5*1e-3];
        hca.XGrid = 'on';
        hca.YGrid = 'on';
        legend(hca,{'entire diff','included for correlation'})
        irf_legend(hca,{['\Sigma_v' sprintf('|<f^{mod}>-<f^{fpi}>|^2 = %g',corr_f_tmp)]},[0.02 0.98])
        pause
      end
    end
  end
end

%% Plot results of correlation analysis
% Get from bgk_timeseries:
% vph_all = [-9000*1e3 -10000*1e3];
% phi_mult_all = [1.0 1.5];
% iff_all = 1:17;
figure(43)
% Ordering: n_vph_all, n_phi_mult_all, n_iff_all
[VPH_ALL,PHI_MULT_ALL,IFF_ALL] = ndgrid(vph_all,phi_mult_all,iff_all);
n_cominations = n_vph_all*n_phi_mult_all*n_iff_all;
n_cominations_vph_phi = n_vph_all*n_phi_mult_all;
h = setup_subplots(2,2,1);
isub = 1;

i_iff_plot = [1 3:n_iff_all]; % remove outlier
colors = pic_colors('matlab');

if 1
  hca = h(isub); isub = isub + 1;
  i_count = 0;
  for i_vph = 1:n_vph_all
    for i_phi_mult = 1:n_phi_mult_all
      i_iff = i_iff_plot;
      i_count = i_count + 1;
      if i_count == 1, hold(hca,'off');
      elseif i_count == 2, hold(hca,'on'); end
      plot_x = corr_n_all(i_vph,i_phi_mult,i_iff); 
      plot_y = corr_j_all(i_vph,i_phi_mult,i_iff);      
      %plot_s = IFF_ALL(ip_vph,i_phi_mult,:);
      %plot_c = abs(VPH_ALL(ip_vph,i_phi_mult,:));
      plot_color = colors(i_count,:);
      plot(hca,plot_x(:),plot_y(:),'.','MarkerSize',15,'color',plot_color,'DisplayName',sprintf('v_{ph} = %.0f km/s, phi _{multiplier} = %.0f',vph_all(i_vph)*1e-3,phi_mult_all(i_phi_mult)))      
      
      hca.XLabel.String = '|n^{mod}-n^{obs}|^2';
      hca.YLabel.String = '|j^{mod}-j^{edi}|^2';
  %pause 
      %scatter(corr_n_all(:),corr_j_all(:),IFF_ALL(:),abs(VPH_ALL(:)))  
    end
  end
  hold(hca,'off')
  legend(hca,'location','northwest')
end
if 1
  hold_on = 0;
  hca = h(isub); isub = isub + 1;
  i_count = 0;
  h_to_leg = hca; h_to_leg(1) = []; % make empty axes array
  for i_vph = 1:n_vph_all
    for i_phi_mult = 1:n_phi_mult_all
      i_count = i_count + 1;    
      for i_iff = i_iff_plot %1:n_iff_all  
        if hold_on == 0, hold(hca,'on'); hold_on = 1; end
        plot_x = corr_n_all(i_vph,i_phi_mult,i_iff); 
        plot_y = corr_j_all(i_vph,i_phi_mult,i_iff);      
        %plot_s = IFF_ALL(ip_vph,i_phi_mult,:);
        %plot_c = abs(VPH_ALL(ip_vph,i_phi_mult,:));        
        plot_color = colors(i_count,:);
        hplot = plot(hca,plot_x(:),plot_y(:),'.','MarkerSize',15,'color',plot_color,'DisplayName',sprintf('v_{ph} = %.0f km/s, phi _{multiplier} = %.0f',vph_all(i_vph)*1e-3,phi_mult_all(i_phi_mult)));
        htext = text(hca,plot_x,plot_y,sprintf('%.0f',i_iff));
        htext.VerticalAlignment = 'top';
        %h_to_leg(i_vph,i_phi_mult) = hplot;
        
        hca.XLabel.String = '|n^{mod}-n^{obs}|^2';
        hca.YLabel.String = '|j^{mod}-j^{edi}|^2';
      end
    end
  end
  hca.Box = 'on';
  %hca.XGrid = 'on';
  %hca.YGrid = 'on';
  hold(hca,'off')
  %legend(h_to_leg,'location','eastoutside')
end
if 0
  hca = h(isub); isub = isub + 1; hold(hca,'off')
  plot_x = VPH_ALL(:,1,:); 
  plot_y = IFF_ALL(:,1,:);
  plot_s = corr_j_all(:,1,:);
  plot_c = plot_s;
  scatter(plot_x(:),plot_y(:),plot_s(:),plot_c(:))
end
if 1
  hca = h(isub); isub = isub + 1;
  i_count = 0;
  for i_vph = 1:n_vph_all
    for i_phi_mult = 1:n_phi_mult_all 
      i_iff = i_iff_plot;           
      i_count = i_count + 1;
      if i_count == 1, hold(hca,'off');
      elseif i_count == 2, hold(hca,'on'); end
      plot_x = corr_f_all(i_vph,i_phi_mult,i_iff); 
      plot_y = corr_j_all(i_vph,i_phi_mult,i_iff);      
      %plot_s = IFF_ALL(ip_vph,i_phi_mult,:);
      %plot_c = abs(VPH_ALL(ip_vph,i_phi_mult,:));
      plot_color = colors(i_count,:);
      plot(hca,plot_x(:),plot_y(:),'.','MarkerSize',15,'color',plot_color,'DisplayName',sprintf('v_{ph} = %.0f km/s, phi _{multiplier} = %.0f',vph_all(i_vph)*1e-3,phi_mult_all(i_phi_mult)))      
      
      hca.XLabel.String = '|f^{mod}-f^{obs}|^2';
      hca.YLabel.String = '|j^{mod}-j^{edi}|^2';
  %pause 
      %scatter(corr_n_all(:),corr_j_all(:),IFF_ALL(:),abs(VPH_ALL(:)))  
    end
  end
  hold(hca,'off')
  legend(hca,'location','northwest')
end
if 1
  hold_on = 0;
  hca = h(isub); isub = isub + 1;
  i_count = 0;
  h_to_leg = hca; h_to_leg(1) = []; % make empty axes array
  for i_vph = 1:n_vph_all
    for i_phi_mult = 1:n_phi_mult_all
      i_count = i_count + 1;    
      for i_iff = i_iff_plot %1:n_iff_all  
        if hold_on == 0, hold(hca,'on'); hold_on = 1; end
        plot_x = corr_f_all(i_vph,i_phi_mult,i_iff); 
        plot_y = corr_j_all(i_vph,i_phi_mult,i_iff);      
        %plot_s = IFF_ALL(ip_vph,i_phi_mult,:);
        %plot_c = abs(VPH_ALL(ip_vph,i_phi_mult,:));        
        plot_color = colors(i_count,:);
        hplot = plot(hca,plot_x(:),plot_y(:),'.','MarkerSize',15,'color',plot_color,'DisplayName',sprintf('v_{ph} = %.0f km/s, phi _{multiplier} = %.0f',vph_all(i_vph)*1e-3,phi_mult_all(i_phi_mult)));
        htext = text(hca,plot_x,plot_y,sprintf('%.0f',i_iff));
        htext.VerticalAlignment = 'top';
        %h_to_leg(i_vph,i_phi_mult) = hplot;
        
        hca.XLabel.String = '|f^{mod}-f^{obs}|^2';
        hca.YLabel.String = '|j^{mod}-j^{edi}|^2';
      end
    end
  end
  hca.Box = 'on';
  %hca.XGrid = 'on';
  %hca.YGrid = 'on';
  hold(hca,'off')
  %legend(h_to_leg,'location','eastoutside')
end
%plot(phi_mult_all)


%%


