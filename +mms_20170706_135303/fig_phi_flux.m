
% Get potential from Efield
tint = irf.tint('2017-07-06T13:54:05.490Z/2017-07-06T13:54:05.620Z');
vph = -9000*1e3;        

c_eval('Epar = gseE?par;',mms_id)
[phi,phi_progressive,phi_ancillary] = get_phi(Epar,vph,tint,tint);

colors = mms_colors('matlab'); 
h = irf_plot(3);

hca = irf_panel('phi');
irf_plot(hca,phi);
hca.YLabel.String = {'\phi','(V)'};
hca.YLabel.Interpreter = 'tex';
irf_legend(hca,sprintf('v_{ph}=%.0f km/s',vph*1e-3),[0.02 0.98],'k')

hca = irf_panel('flux edi');
flux_scale = 1e-6;
irf_plot(hca,{flux0_mms1.tlim(tint)*flux_scale,flux180_mms1.tlim(tint)*flux_scale},'comp');
hca.YLabel.String = {'flux','(10^{6} cm^{-1}s^{-1}sr^{-1})'};
hca.YLabel.Interpreter = 'tex';
irf_legend(hca,{'0^{o}','180^{o}'},[0.02 0.98])

hca = irf_panel('phi, flux edi 180');
irf_plot(hca,phi*-1);
set(hca,'YColor',colors(1,:)); % color of axis lines and numbers
hca.YLabel.String = {'\phi','(V)'};

ax2 = axes('Position',get(hca,'Position'));
irf_plot(ax2,flux180_mms1.tlim(tint)*flux_scale,'color',colors(2,:));
set(ax2,'XAxisLocation','top','xtick',[],'xlabel',[]); % remove 'xtick' if xticks required
set(ax2,'YAxisLocation','right');
set(ax2,'Color','none'); % color of axis
set(ax2,'YColor',colors(2,:)); % color of axis lines and numbers


irf_zoom(h,'x',tint)
irf_zoom(ax2,'x',tint)

hca.YLabel.String = {'-\phi','(V)'};
hca.YLabel.Interpreter = 'tex';

ax2.YLabel.String = {'flux','(10^{6} cm^{-1}s^{-1}sr^{-1})'};
ax2.YLabel.Interpreter = 'tex';
set(ax2,'XAxisLocation','top','xtick',[],'xlabel',[]); % remove 'xtick' if xticks required

h(1).XGrid = 'off';
h(2).XGrid = 'off';
h(3).XGrid = 'off';
ax2.XGrid = 'off';
h(1).YGrid = 'off';
h(2).YGrid = 'off';
h(3).YGrid = 'off';
ax2.YGrid = 'off';


