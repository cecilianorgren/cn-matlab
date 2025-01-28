
timeOverview = irf.tint('2015-10-16T10:33:40.00Z/2015-10-16T10:33:50.00Z'); % knee
%tint = timeOverview + [+5 -5];
tintZoom = irf.tint('2015-10-16T10:33:40.00Z/2015-10-16T10:33:50.00Z'); % knee
tint = tintZoom;
%%
wavE1 = irf_wavelet(dslE1brst.abs.tlim(tint),'wavelet_width',5.36,'f',[1 4000]);
wavE1.f_units = 'Hz'; wavE1.f_label = 'f [Hz]';

%%
hca=subplot(2,1,2);
irf_spectrogram(hca,wavE1_par_w2,'log')
hca.YScale = 'log';
hold(hca,'on')
hh=irf_plot(hca,{fce1brst,fcp1brst,flh1brst,fpp1brst,fpe1brst},'comp');
legend(hh.Children(1:5),{'f_{ce}','f_{ci}','f_{LH}','f_{pi}','f_{pe}'},'location','bestoutside');

%%
hca=subplot(2,1,2);
irf_spectrogram(hca,wavE1x_w4,'log')
hca.YScale = 'log';
hold(hca,'on')
hh=irf_plot(hca,{fce1brst,fcp1brst,flh1brst,fpp1brst,fpe1brst},'comp');
legend(hh.Children(1:5),{'f_{ce}','f_{ci}','f_{LH}','f_{pi}','f_{pe}'},'location','bestoutside');

%%
hca=subplot(2,1,1);
cb(1)=colorbar
hca.CLim = [-14 -2]
hca=subplot(2,1,2);
cb(2)=colorbar
hca.CLim = [-14 -2]

%irf_plot(hca,dmpaB1brst)

%%
ic = 1;
fn = 100;
c_eval('wavE?x = irf_wavelet(dslE?brst.x.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''fn'',100);',ic)
c_eval('wavE?x.f_units = ''Hz''; wavE?x.f_label = ''f [Hz]''; wavE?x.p_label = {''log_{10} E_x^2'',''(mV/m)^2/Hz''};',ic)


c_eval('wavE?y = irf_wavelet(dslE?brst.y.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''fn'',100);',ic)
c_eval('wavE?y.f_units = ''Hz''; wavE?y.f_label = ''f [Hz]''; wavE?y.p_label = {''log_{10} E_y^2'',''(mV/m)^2/Hz''};',ic)

c_eval('wavB?x = irf_wavelet(Bscm.x.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''fn'',100);',ic)
c_eval('wavB?x.f_units = ''Hz''; wavB?x.f_label = ''f [Hz]''; wavB?x.p_label = {''log_{10} B_x^2'',''(nT)^2/Hz''};',ic)

%%
ic = 1:4;
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:38.00Z')+[-2 2];
%c_eval('dtB = dmpaB?scm.time(2)-dmpaB?scm.time(1); newTime = dmpaB?scm.time(1):dtB/4:dmpaB?scm.time(end); wavB? = irf_wavelet(dmpaB?scm.abs.resample(newTime).tlim(tint),''wavelet_width'',5.36*4,''f'',[1 4000],''nf'',100);',ic)
%c_eval('dtE = dslE?brst.time(2)-dslE?brst.time(1); newTime = dslE?brst.time(1):dtE/4:dslE?brst.time(end); wavE? = irf_wavelet(dslE?brst.abs.resample(newTime).tlim(tint),''wavelet_width'',5.36*4,''f'',[1 4000],''nf'',100);',ic)
c_eval('wavE? = irf_wavelet(dslE?brst.abs.tlim(tint),''wavelet_width'',5.36,''f'',[1 5000],''nf'',100);',ic)
c_eval('wavB? = irf_wavelet(dmpaB?scm.abs.tlim(tint),''wavelet_width'',5.36,''f'',[1 5000],''nf'',100);',ic)

c_eval('wavB?.p =wavB?.p{1};',ic); 
c_eval('wavE?.p =wavE?.p{1};',ic); 
c_eval('wavE?.f_units = ''Hz''; wavE?.f_label = ''f [Hz]''; wavE?.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)
c_eval('wavB?.f_units = ''Hz''; wavB?.f_label = ''f [Hz]''; wavB?.p_label = {''log_{10} B^2'',''(nT)^2/Hz''};',ic)
%% 
h = irf_plot(8);
ic = 1;
hca = irf_panel('B');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {'B','(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('brst E');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y},''comp'');',ic)
hca.YLabel.String = {'E','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_x','E_y'},[0.95 0.95]);
irf_zoom(hca,'y')

hca=irf_panel('edist');
c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E_{e,OMNI}','(eV)'},'Interpreter','tex');

if 0
hca=irf_panel('edistper');
c_eval('irf_spectrogram(hca,ePSDperp?,''log'',''donotfitcolorbarlabel'');',ic)
%irf_legend(hca,'(f)',[0.99 0.98],'color','k','fontsize',12)
set(hca,'yscale','log');
hca.CLim = [-2 5];
%caxis(hca,1*[-1 1])
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E','(eV)'});
%colormap(hca,cn.cmap('bluered3'))
end

hca=irf_panel('Ex wavelet');
c_eval('irf_spectrogram(hca,wavE?x,''log'',''donotfitcolorbarlabel'');',ic)
%hcb(1) = colorbar('peer',hca);
%hcb(1).YLabel.String = 'E_x [(mV/m)^2/Hz]';
hca.YScale = 'log';
hold(hca,'on')
c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
%legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
hca.YLabel.String = 'f [Hz]';
hca.CLim = [-9 -2.5]; 
hca.YTick = [1e2 1e3 1e4];

hca=irf_panel('Ey wavelet');
c_eval('irf_spectrogram(hca,wavE?y,''log'',''donotfitcolorbarlabel'');',ic)
%hcb(2) = colorbar('peer',hca);
%hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
hca.YScale = 'log';
hold(hca,'on')
c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
%legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
hca.YLabel.String = 'f [Hz]';
hca.CLim = [-9 -2.5]; 
hca.YTick = [1e2 1e3 1e4];

hca=irf_panel('Bx wavelet');
c_eval('irf_spectrogram(hca,wavB?,''log'',''donotfitcolorbarlabel'');',ic)
%hcb(2) = colorbar('peer',hca);
%hcb(2).YLabel.String = 'E_y [(mV/m)^2/Hz]';
hca.YScale = 'log';
hold(hca,'on')
c_eval('hfce=irf_plot(hca,fce?brst,''cyan'');',ic)
c_eval('hfpp=irf_plot(hca,fpp?brst,''green'');',ic)
irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
%legend(hca,[hfce hfpp],{'f_{ce}','f_{pi}'},'location','northwest');
hca.YLabel.String = 'f [Hz]';
hca.CLim = [-9 -2.5]; 
hca.YTick = [1e2 1e3 1e4];

clim = [-1 1];
ylim = [10 1e3];

if 1
  
    hca = irf_panel(irf_ssub('edistparperp ?',ic));
    c_eval('[~,hcb] = irf_spectrogram(hca,ePSDparperp?,''log'',''donotfitcolorbarlabel'');',ic)
    hold(hca,'on')  
    hold(hca,'off')
    %irf_legend(hca,'(g)',[0.99 0.98],'color','k','fontsize',12)
    set(hca,'yscale','log');
    hca.CLim = clim;
    set(hca,'ytick',[1e1 1e2 1e3 1e4]);
    hca.YLim = ylim;
    ylabel(hca,{'E_e','(eV)'});
    colormap(hca,cn.cmap('bluered3'))
    
    hhcb(isub) = hcb;
    if 1
      hcb.YTick = 0.6*[-1 1];
      hcb.YTickLabel = {'f_{perp}','f_{a/par}'};
      hcb.YLabel.String = irf_ssub('mms {?}',ic);
      hcb.YLabel.Color = mms_colors(irf_ssub('?',ic));
    end
end
if 1
  hca = irf_panel(irf_ssub('edistparapar ?',ic));
  c_eval('[~,hcb] = irf_spectrogram(hca,ePSDparapar?,''log'',''donotfitcolorbarlabel'');',ic)
  %irf_legend(hca,'(f)',[0.99 0.98],'color','k','fontsize',12)
  set(hca,'yscale','log');
  hca.CLim = clim;
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLim = ylim;
  ylabel(hca,{'E','(eV)'});
  colormap(hca,cn.cmap('bluered3'))

  %hca.YLabel.String = {irf_ssub('E_{?}',ic),'(eV)'};

  % Colorbar labels
  if 1
    hcb.YTick = 0.6*[-1 1];
    hcb.YTickLabel = {'f_{apar}','f_{par}'};
    hcb.YLabel.String = irf_ssub('mms {?}',ic);
    hcb.YLabel.Color = mms_colors(irf_ssub('?',ic));
  end
end

irf_plot_axis_align
irf_zoom(h,'x',tint)
h(1).Title.String = irf_ssub('MMS ?',ic)