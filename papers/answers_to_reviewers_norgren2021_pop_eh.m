h = irf_plot(6);
ffilt = 100;

hca = irf_panel('Epar');
irf_plot(hca,gseE1par)
hca.YLabel.String = 'E_{||} (mV/m)';
hca.YLabel.Interpreter = 'tex';

hca = irf_panel('phi');
irf_plot(hca,phi1)
hca.YLabel.String = '\phi (V)';
hca.YLabel.Interpreter = 'tex';

hca = irf_panel('scpot filt');
%ffilt = 100;
irf_plot(hca,scPot1.filt(ffilt,0,[],5))
hca.YLabel.String = '\delta V_{sc} (V)';
hca.YLabel.Interpreter = 'tex';
irf_legend(hca,sprintf('f>%g Hz',ffilt),[0.05 0.98])

hca = irf_panel('scpot');
irf_plot(hca,scPot1)
hca.YLabel.String = 'V_{sc} (V)';
hca.YLabel.Interpreter = 'tex';


hca = irf_panel('scpot filt / scpot');
irf_plot(hca,scPot1.filt(ffilt,0,[],5)/scPot1)
hca.YLabel.String = '\delta V_{sc}/V_{sc}';
hca.YLabel.Interpreter = 'tex';

if 0
  hca = irf_panel('n scpot');
  %ffilt = 100;
  irf_plot(hca,nescpot.filt(ffilt,0,[],5)*1e3)
  hca.YLabel.String = 'n_{V_{sc}} (10^{-3} cm^{-3})';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,sprintf('f>%g Hz',ffilt),[0.1 0.98])
end
if 0
  hca = irf_panel('n scpot norm');
  %ffilt = 100;
  irf_plot(hca,nescpot.filt(ffilt,0,[],5)/nescpot)
  hca.YLabel.String = '\delta V_{sc}/V_{sc}';
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('f>%g Hz',ffilt),[0.1 0.98])
end

hca = irf_panel('dcv');
dcv1_spinplane = dcv1.clone(dcv1.time,dcv1.data(:,1:4)); 
irf_plot(hca,dcv1_spinplane)
hca.YLabel.String = 'V_{p-sc} (V)';
hca.YLabel.Interpreter = 'tex';
irf_legend(hca,'spin plane probes',[0.02 0.05])
%c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
legends_color = {'k','k','k','k','k','k','k','k','k','k','k','k'};
for ipanel = 1:numel(h)
  irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
  h(ipanel).YLabel.FontSize = fontsize;
  h(ipanel).FontSize = fontsize;
end

irf_plot_axis_align
irf_zoom(h,'x',phi1.time); 
irf_zoom(h,'y')