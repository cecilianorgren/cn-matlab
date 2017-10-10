% Load electron distributions


tint = irf.tint('2015-10-16T10:33:15.00Z/2015-10-16T10:33:50.00Z');

% Load B field
c_eval('Bvec? = dmpaB1/dmpaB1.abs;',ic);
c_eval('Bvec = Bvec?;',ic);

% ExB velocity in eV
mp = 1.6726e-27;
e = 1.6022e-19;
c_eval('vExBeV? = vExB?.abs*vExB?.abs*1e6*0.5*mp/e;',ic);
%c_eval('vExBeV? = mms.fftbandpass(vExBeV?,0,0.5);',ic);

% Ion velocity in eV
c_eval('vieV? = vi?brst.abs*vi?brst.abs*1e6*0.5*mp/e;',ic);

% Electron velocity in eV
c_eval('veeV? = ve?brst.abs*ve?brst.abs*1e6*0.5*mp/e;',ic);


h=irf_plot(9,'newfigure');

c_eval('dmpaB = dmpaB?.tlim(tint);',ic)
hca=irf_panel('B');
irf_plot(hca,dmpaB);
ylabel(hca,'B_{DMPA} (nT)','Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.5 0.1])
irf_legend(hca,'(a)',[0.99 0.98],'color','k')

hca=irf_panel('ve');
c_eval('irf_plot(hca,ve?brst);',ic)
ylabel(hca,'v_{e} (km s^{-1})','Interpreter','tex');
irf_legend(hca,{'V_{x}','V_{y}','V_{z}'},[0.5 0.1])
irf_legend(hca,'(b)',[0.99 0.98],'color','k')

hca = irf_panel('brst E');
c_eval('irf_plot(hca,{dslE?brst.x,dslE?brst.y},''comp'');',ic)
hca.YLabel.String = 'E_{DSL} [mV/m]';
irf_legend(hca,{'E_x','E_y'},[0.95 0.95]);

hca = irf_panel('brst scPot');
c_eval('irf_plot(hca,(-1)*P?brst);',ic);
hca.YLabel.String = '-scPot [V]';


hca=irf_panel('n');
c_eval('irf_plot(hca,ne?brst);',ic)
hold(hca,'on');
c_eval('irf_plot(hca,ni?brst,''b'');',ic)
hold(h(3),'off');
set(hca,'yscale','lin');
irf_zoom(hca,'y',[1 300]);
ylabel(hca,'n (cm^{-3})','Interpreter','tex');
irf_legend(hca,{'n_{e}','n_{i}'},[0.5 0.9]);
irf_legend(hca,'(c)',[0.99 0.98],'color','k')

hca=irf_panel('idist');
irf_spectrogram(hca,specfiomni,'log','donotfitcolorbarlabel');
hold(hca,'on');
c_eval('irf_plot(hca,vieV?);',ic);
c_eval('irf_plot(hca,vExBeV?,''w'');',ic);
hold(hca,'off');
irf_legend(hca,'(d)',[0.99 0.98],'color','w','fontsize',12);
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
irf_legend(hca,{'v_i'},[0.5 0.9]);
irf_legend(hca,'v_{ExB}',[0.6 0.9],'color','w');
ylabel(hca,'E_{i} (eV)','fontsize',12,'Interpreter','tex');

hca=irf_panel('edist');
irf_spectrogram(hca,specfomni,'log','donotfitcolorbarlabel');
hold(hca,'on');
c_eval('irf_plot(hca,Te?brst)',ic)
hold(hca,'off');
irf_legend(hca,'(e)',[0.99 0.98],'color','w','fontsize',12)
irf_legend(hca,{'T_{x}','T_{y}','T_{z}'},[0.5 0.1])
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
%caxis(hca,[-15 -5])
ylabel(hca,'E_{e} (eV)','fontsize',12,'Interpreter','tex');

hca=irf_panel('edistparperp');
irf_spectrogram(hca,specparperp,'log','donotfitcolorbarlabel');
irf_legend(hca,'(g)',[0.99 0.98],'color','k','fontsize',12)
set(hca,'yscale','log');
caxis(hca,[-2 2])
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,'E (eV)','fontsize',12);
colormap(hca,cn.cmap('bluered3'))

hca=irf_panel('edistparapar');
irf_spectrogram(hca,specparapar,'log','donotfitcolorbarlabel');
irf_legend(hca,'(f)',[0.99 0.98],'color','k','fontsize',12)
set(hca,'yscale','log');
caxis(hca,[-2 2])
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,'E (eV)','fontsize',12);
colormap(hca,cn.cmap('bluered3'))


title(h(1),'MMS 1')

%tints = irf.tint('2015-10-01T06:53:30.00Z/2015-10-01T06:54:10.00Z');
irf_plot_axis_align(h(1:7));
irf_zoom(h(1:7),'x',tint);
set(h(1:7),'fontsize',12);

for ii = 1: numel(h);
    h(ii).Position(3) = h(ii).Position(3)*0.85;
    grid(h(ii),'off')
end

%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
%print('-dpng','-painters','-r600','overviewedists1s.png');