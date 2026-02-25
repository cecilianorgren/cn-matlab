h = setup_subplots(1,2);
isub = 1;
hca = h(isub); isub = isub + 1;
vdf = pdist.reduce('2D',[1 0 0],[0 0 1]);
vdf.plot_plane(hca);
hca.XLabel.String = 'v_x (km/s)';
hca.YLabel.String = 'v_z (km/s)';

hca = h(isub); isub = isub + 1;
%scatter(hca,MP.vx,MP.vz,10,MP.dv.*MP.df,'.')
scatter(hca,MP.vx,MP.vz,2,'o','filled')
hca.Box = 'on';
hca.XLabel.String = 'v_x (km/s)';
hca.YLabel.String = 'v_z (km/s)';

c_eval('axis(h(?),''square'')',1:numel(h))
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).FontSize = 12;',1:numel(h))
compact_panels(h,0.09,0.25,1)

hlinks = linkprop(h,{'XLim','YLim'});
hca.XLim = 3000*[-1 1];
hca.YLim = 3000*[-1 1];