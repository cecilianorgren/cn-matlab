fn=figure(71);h=irf_plot(2);

if 1,   % PANEL: C?       WHISPER spectrogram
    hca=irf_panel('WHISPER spectrogram natural');
    varname=irf_ssub('Electric_Spectral_Power_Density__C?_CP_WHI_NATURAL',ic);
    %[~,~,~,varunits]=c_caa_var_get(varname);
    varunits='(V/m)^2/Hz';
    % REMOVE 'fillspectgrogramgaps' flag in the next line if precise intervals of
    % WHISPER measurements are needed !!!!
    irf_plot(hca,varname,'tint',tint,'colorbarlabel',varunits,'fitcolorbarlabel','fillspectrogramgaps','nolabels');
    hold(hca,'on');
    c_eval('fpe=irf_plasma_calc(irf_resamp(B?,ncal_PEACE?),ncal_PEACE?,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
    irf_plot(hca,irf_tappl(fpe,'/1e3'),'-','linewidth',0.2,'color','k');
    c_eval('fpemom=irf_plasma_calc(irf_resamp(B?,nPEACE),nPEACE,0,0,0,''Fpe'');',ic); % calculate electron gyrofrequency
    irf_plot(hca,irf_tappl(fpemom,'/1e3'),'.','linewidth',0.2,'color','r');
    hold(hca,'off');
    caxis(hca,[-16 -11]);
    set(hca,'yscale','log','ytick',[3 4 5 1e1 20 30 50 ]);
    irf_zoom(hca,'y',[2 12]);
end
if 1
    hca=irf_panel('Electric field')
    irf_plot(hca,gsmE3)
end