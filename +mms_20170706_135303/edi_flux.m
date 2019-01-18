tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z');

figure(51);
mms_id = 1;
npanels = 3;
h = irf_plot(npanels);

if 1
  hca = irf_panel('flux nodes');
  c_eval('irf_plot(hca,flux_135_180_mms?);',mms_id) 
  irf_legend(hca,{'sect. 1','sect. 2','sect. 3','sect. 4'},[0.98 0.98])
  hca.YLabel.String = {'Field-of-view flux',sprintf('(%s)',flux_135_180_mms1.units)};
end
if 1
  hca = irf_panel('flux scaled nodes');
  c_eval('irf_plot(hca,flux_135_180_scaled_mms?);',mms_id) 
  irf_legend(hca,{sprintf('sect. 1 x %.3f',pa_scaling(1)),sprintf('sect. 2 x %.3f',pa_scaling(2)),sprintf('sect. 3 x %.3f',pa_scaling(3)),sprintf('sect. 4 x %.3f',pa_scaling(4))},[0.98 0.98])
  hca.YLabel.String = {'Parallel flux',sprintf('(%s)',flux_135_180_scaled_mms1.units)};
end
if 1
  hca = irf_panel('flux scaled summed');
  c_eval('irf_plot(hca,flux_apar_int_mms?);',mms_id) 
  hca.YLabel.String = {'Total parallel flux',sprintf('(%s)',flux_apar_int_mms1.units)};
end

h(1).Title.String = sprintf('MMS %g',mms_id);
irf_zoom(h,'x',tint_zoom); % from mms_20170706_135303.paper_fig1
irf_zoom(h,'y');

h(1).YLim = [0 8e6];
h(2).YLim = [0 4e6];
h(3).YLim = [0 9e6];

%cn.print(sprintf('integrated_edi_flux_mms%g_sameylims',mms_id))