

h = irf_plot(3);

hca = irf_panel('ve');
irf_plot(hca,mvaVe3)

hca = irf_panel('pressure');
irf_plot(hca,{facPepp3.xx,facPepp3.yy,facPepp3.zz,PB3},'comp')
hleg = irf_legend(hca,{'$p_{e\parallel}$','$p_{e\perp}$','$p_{e\perp}$','$P_B$'},[0.98 0.98]);
c_eval('hleg(?).Interpreter = ''latex'';',1:numel(hleg))

hca = irf_panel('firehose');
irf_plot(hca,{facPepp3.xx-facPepp3.yy,PB3},'comp')
hleg = irf_legend(hca,{'$p_{e\parallel}-p_{e\perp}$','$P_B$'},[0.98 0.98]);
c_eval('hleg(?).Interpreter = ''latex'';',1:numel(hleg))