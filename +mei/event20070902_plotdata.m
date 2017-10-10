% Different plots

%% Overview plot for single SC 
sc=3;
%c_eval('tint=[diB?fgm(1,1) diB?fgm(end,1)];',sc)
tint = tint_ov+[-10 10];
h = irf_plot(7,'newfigure');
isub = 1;

if 1 % Electron differential energy flux
    hca = h(isub); isub=isub+1;
    c_eval('deflux = e_deflux?;',sc)
    
    specrec.t = deflux.t;
    specrec.f = deflux.dep_x{2}.data(1,:)'*1e3; % eV
    specrec.p = {squeeze(nansum(deflux.data,2))};
    irf_spectrogram(hca,specrec.t,specrec.p{1},specrec.f)
    hca.YScale = 'log';
    hca.YTick = [0.01 0.1 1 10 100 1000 10000];
    hca.YLabel.String = 'E_e [eV]';
    hc = colorbar('peer',hca);
    hc.Label.String = 'keV/cm^2 s sr keV';
    
    %irf_plot(hca,irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',sc),'sum_dim1','colorbarlabel','log_{10} dEF\newline #/cm^2 s sr keV','fitcolorbarlabel');
    %hca.YLabel.String = 'E [nT]';
end
if 1 % n (PEACE)
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,irf_ssub('peaNe?',sc));
    hca.YLabel.String = 'n_e [cc]';    
end
if 1 % Te (PEACE)
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,irf_ssub('parTe?',sc));
    hca.YLabel.String = 'T_{e,||} [eV]';    
end
if 1 % beta
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('beta?',sc));
    hca.YLabel.String = '\beta ';
    %irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    %irf_legend(hca,{'GSM'},[0.02, 0.95]);
end

if 0 % B ISR2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('diB?',sc));
    hca.YLabel.String = 'B [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % B GSM
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('gsmB?',sc));
    hca.YLabel.String = 'B [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
try
if 1 % V GSM
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('gsmVi?',sc));
    hca.YLabel.String = 'Vi [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
end
if 0 % E
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('diE?',sc));
    hca.YLabel.String = 'E [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % E
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('gsmE?',sc));
    hca.YLabel.String = 'E [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 0 % E power spectrogram
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftE?.t,fftE?.p{1},fftE?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = '(nT)^2/Hz';
    hca.NextPlot = 'add';
    c_eval('irf_plot(hca,flh?)',sc);
    c_eval('irf_plot(hca,fce?)',sc);
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';
end
if 0 % B power spectrogram
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftB?.t,fftB?.p{1},fftB?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = '(mV/m)^2/Hz';
    hca.NextPlot = 'add';
    irf_plot(hca,irf_ssub('flh?',sc));
    irf_plot(hca,irf_ssub('fce?',sc));
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';    
end
if 0
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftEB?.t,fftEB?.p{1},fftEB?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = 'mV/m/nT/Hz';
    hca.NextPlot = 'add';
    irf_plot(hca,irf_ssub('flh?',sc));
    irf_plot(hca,irf_ssub('fce?',sc));
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';    
end


title(h(1),irf_ssub('C?',sc))
irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align
set(gcf,'paperpositionmode','auto');
strPrint = [datestr(irf_time(tint(1),'epoch>datenum'),'yyyy-mm-dd') '_Overview'];
%eval(['print -depsc /Users/Cecilia/Research/LH2/BurstModeOverview/' strPrint '.eps'])