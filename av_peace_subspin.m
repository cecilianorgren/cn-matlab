%% better than spin resolution Fighighres
if 1, % initialize figure
    fn=figure(61);
    h=irf_plot(6);
end
sclist=3;ic=3;
if 1, % read FGM data form all sc
    %    c_eval('[caaB?,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');');
    c_eval('[caaB?,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_5VPS'');',sclist);
    c_eval('B?=irf_abs(B?);',sclist);
    c_eval('gsmB?=irf_gse2gsm(B?);',sclist);
end
if 1, % read EFW data
    %     c_eval('[caaE?,~,diE?,diE?units]=c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_E'');',sclist);
%     c_eval('[caaE?,~,diE?,diE?units]=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L3_E3D_INERT'');',sclist);
    c_eval('[caaE?,~,diE?,diE?units]=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'');',sclist);
end
if 0,   % PANEL: FGM B single s/c
    hca=irf_panel('B');
    c_eval('irf_plot(hca,gsmB?);',ic);
    ylabel(hca,'B [nT] GSM');
    irf_zoom(hca,'y','smart');
    irf_legend(hca,{'B_X','B_Y','B_Z','B'},[0.02 0.3])
    irf_legend(hca,{['C' num2str(ic)]},[0.98 0.98],'color','k')
end
if 1,   % PANEL: EFW E field in ISR2 reference frame single s/c
    hca=irf_panel('E');
    c_eval('irf_plot(hca,diE?)',ic);
    ylabel(hca,'E [mV/m] ISR2');
    irf_zoom(hca,'ylim','smart');
    irf_legend(hca,{'E_X','E_Y'},[0.02 0.19])
    irf_legend(hca,{['C' num2str(ic)]},[0.98 0.98],'color','k')
end
if 0,   % PANEL: CIS HIA single s/c spectrogam
    hca=irf_panel('CIS');
    if ic~=2,
        irf_plot(hca,['flux__C' num2str(ic) '_CP_CIS_HIA_HS_1D_PEF'],'colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis(hca,[3.9 6.1]);
        set(hca,'yscale','log')
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel(hca,'E [eV]')
        irf_legend(hca,{['C' num2str(ic)]},[0.98 0.05],'color','k')
        irf_colormap('default');
    end
end
if 0,   % PANEL: PEACE single s/c spectrogram
    hca=irf_panel('PEACE spec');
    irf_plot(hca,['Data__C' num2str(ic) '_CP_PEA_PITCH_SPIN_DPFlux'],'sum_dim1','colorbarlabel','log10 dPF\newline 1/cm^2 s sr keV','fitcolorbarlabel');
    caxis(hca,[5.9 7.6]);
    set(hca,'yscale','log','ylim',[100 3e4])
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    irf_legend(hca,{['C' num2str(ic)]},[0.98 0.05],'color','k')
    ylabel('E [eV]');
end
if 0,   % PANEL: PEACE PEA_PITCH_3DRL_PSD high res C3
    hca=irf_panel('PEACE PITCH 3DRL');
    ic=3;
    res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_PITCH_3DRL_PSD',ic));
    [delmett,ind]=irf_tlim(res.tt,tint);
    specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
    if 0, % energy spectorgram (integrated over pitch angles)
        specrec.f=log10(res.en);
        specrec.p=res.omni(ind,:);
        specrec.f_label=['Log10 ' res.enlabel];
    elseif 1, % pitch angle spectrogram for given energy
        specrec.f=res.theta;specrec.f_label='Pitch angle';
        specrec.p=res.pitch_angle(ind,:);
        enindex=13;
        res.en(enindex)
        specrec.f_label=[specrec.f_label '  \newline[E=' num2str(res.en(enindex),4) 'eV]'];
        specrec.p=log10(res.data(ind,:,enindex));
    end
    specrec_C33DRH=specrec;
    irf_spectrogram(hca,specrec);
    caxis(hca,[-1.99 0.49])
    irf_legend(hca,['C' num2str(ic)],[0.98,0.98]);
    set(hca,'ytick',[30 60 90 120 150]);
end
if 1,   % PANEL: PEACE 3DXPH_DEFlux high res energy spectrogram
    hca=irf_panel('PEACE 3DXPH_DEFlux energy');
    res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_3DXPH_PSD',ic));
    [~,ind]=irf_tlim(res.tt,tint);
    specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
    if 1, % energy spectorgram (integrated over pitch angles)
        specrec.f=log10(res.en);
        specrec.p=res.omni(ind,:);
        specrec.f_label=['Log10 ' res.enlabel];
        irf_spectrogram(hca,specrec);
    elseif 0, % pitch angle spectrogram for given energy
        specrec.f=res.theta;specrec.f_label='Pitch angle';
        specrec.p=res.pitch_angle(ind,:);
        enindex=13;
        specrec.f_label=[specrec.f_label '  \newline[E=' num2str(res.en(enindex),4) 'eV]'];
        specrec.p=log10(res.data(ind,:,enindex));
        irf_spectrogram(hca,specrec);
        set(hca,'ytick',[30 60 90 120 150]);
    end
    caxis(hca,[-1.99 0.49]);
    irf_legend(hca,['C' num2str(ic)],[0.98,0.98]);
end
if 1,   % PANEL: PEACE 3DXPH_DEFlux high res angular spectrogra,
    hca=irf_panel('PEACE 3DXPH_DEFlux angular');
    res=c_caa_construct_subspin_res_data(irf_ssub('Data__C?_CP_PEA_3DXPH_PSD',ic));
    [delmett,ind]=irf_tlim(res.tt,tint);
    specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PSD [' res.dataunits ']']);
    if 0, % energy spectorgram (integrated over pitch angles)
        specrec.f=log10(res.en);
        specrec.p=res.omni(ind,:);
        specrec.f_label=['Log10 ' res.enlabel];
        irf_spectrogram(hca,specrec);
    elseif 1, % pitch angle spectrogram for given energy
        specrec.f=res.theta;specrec.f_label='Pitch angle';
        specrec.p=res.pitch_angle(ind,:);
        enindex=14; % specify which energy chanel
        specrec.f_label=[specrec.f_label '  \newline[E=' num2str(res.en(enindex),4) 'eV]'];
        specrec.p=log10(res.data(ind,:,enindex));
        irf_spectrogram(hca,specrec);
        set(hca,'ytick',[30 60 90 120 150]);
    end
    caxis(hca,[-1.99 0.49]);
    irf_legend(hca,['C' num2str(ic)],[0.98,0.98]);
end
if 1,  % PANEL: CIS HIA high res energy C3
    hca=irf_panel('CIS CODIF high res energy');
    ic=3;
    %res=c_caa_construct_subspin_res_data(irf_ssub('x3d_ions__C?_CP_CIS_CODIF_HS_H1_PSD',ic));
    res=c_caa_construct_subspin_res_data(irf_ssub('x3d_ions__C?_CP_CIS_HIA_HS_MAG_IONS_PEF',ic));
    
    [~,ind]=irf_tlim(res.tt,tint);
    specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PEF \newline [' res.dataunits ']']);
    if 1, % energy spectorgram (integrated over pitch angles)
        specrec.f=log10(res.en);
        specrec.p=res.omni(ind,:);
        specrec.f_label='Log_{10} E [eV]';
        yticks=[1 2 3 4 5];
    elseif 1, % pitch angle spectrogram for given energy
        specrec.f=res.theta;specrec.f_label='Pitch angle';
        specrec.p=res.pitch_angle(ind,:);
        enindex=(26:30);
        if numel(enindex)==1,
            specrec.f_label=[specrec.f_label '\newline[E=' num2str(res.en(enindex),'%6.f') 'eV]'];
        else
            specrec.f_label=[specrec.f_label '\newline[E=' num2str(res.en(enindex(1)),'%6.f') ' - ' num2str(res.en(enindex(end)),'%6.f') 'eV]'];
        end
        specrec.p=sum(res.data(ind,:,enindex),3);
        yticks=[45 90 135 ];
    end
    irf_spectrogram(hca,specrec);
    caxis([1 5]);
    set(hca,'ytick',yticks);
    irf_legend(hca,['C' num2str(ic)],[0.02,0.98]);
end
if 0,  % PANEL: CIS CODIF high res energy C4
    hca=irf_panel('CIS CODIF high res pitch angle');
    %res=c_caa_construct_subspin_res_data(irf_ssub('x3d_ions__C?_CP_CIS_CODIF_HS_H1_PSD',ic));
    res=c_caa_construct_subspin_res_data(irf_ssub('x3d_ions__C?_CP_CIS_HIA_HS_MAG_IONS_PEF',ic));
    
    [~,ind]=irf_tlim(res.tt,tint);
    specrec=struct('t',res.tt(ind),'dt',res.dtsampling/2,'p_label',['Log PEF  \newline [' res.dataunits ']']);
    if 0, % energy spectorgram (integrated over pitch angles)
        specrec.f=log10(res.en);
        specrec.p=res.omni(ind,:);
        specrec.f_label='Log_{10} E [eV]';
        yticks=[1 2 3 4 5];
    elseif 1, % pitch angle spectrogram for given energy
        specrec.f=res.theta;specrec.f_label='Pitch angle';
        specrec.p=res.pitch_angle(ind,:);
        enindex=(26:30);
        if numel(enindex)==1,
            specrec.f_label=[specrec.f_label '\newline[E=' num2str(res.en(enindex),'%6.f') 'eV]'];
        else
            specrec.f_label=[specrec.f_label '\newline[E=' num2str(res.en(enindex(1)),'%6.f') ' - ' num2str(res.en(enindex(end)),'%6.f') 'eV]'];
        end
        specrec.p=sum(res.data(ind,:,enindex),3);
        yticks=[45 90 135 ];
    end
    irf_spectrogram(hca,specrec);
    caxis([4 7]);
    set(hca,'ytick',yticks);
    irf_legend(hca,['C' num2str(ic)],[0.02,0.98]);
end
irf_plot_axis_align
irf_zoom(h,'x',tint);
irf_legend(h(1),'Fighighres',[0 1.001],'fontsize',8,'color',[0.5 0.5 0.5]);    irf_plot_axis_align
irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);
irf_timeaxis(h);
