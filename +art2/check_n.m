%% load data
cd /Users/Cecilia/Data/BM/20070831/
tint = toepoch([2007 08 31 10 13 00;2007 08 31 10 25 00])';
% B-field
c_eval('gseB?=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3:4);
c_eval('diB?=c_coord_trans(''GSE'',''DSI'',gseB?,''cl_id'',?);',3:4);

% E-field
c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',3:4);

% E-field spectra
c_eval('specE?.p=c_caa_var_get(''ESUM__C?_CP_STA_PPP'',''mat'');',3:4);
c_eval('specE?.t=c_caa_var_get(''Time__C?_CP_STA_PPP'',''mat'');',3:4);
c_eval('specE?.f=c_caa_var_get(''Frequency__C?_CP_STA_PPP'',''mat'');',3:4);
c_eval('specE?.p=specE?.p(:,[2:end]);',3:4) % take away time row
% irf_spectrogram(specE3)

c_eval('specEpol?.p=c_caa_var_get(''POLSVD__C?_CP_STA_PPP'',''mat'');',3:4);
c_eval('specEpol?.t=c_caa_var_get(''Time__C?_CP_STA_PPP'',''mat'');',3:4);
c_eval('specEpol?.f=c_caa_var_get(''Frequency__C?_CP_STA_PPP'',''mat'');',3:4);
c_eval('specEpol?.p=specEpol?.p(:,[2:end]);',3:4) % take away time row
%c_eval('specEpol?.p(specEpol?.p<0) = NaN;',3:4)

c_eval('specEell?.p=c_caa_var_get(''ELLSVD__C?_CP_STA_PPP'',''mat'');',3:4);
c_eval('specEell?.t=c_caa_var_get(''Time__C?_CP_STA_PPP'',''mat'');',3:4);
c_eval('specEell?.f=c_caa_var_get(''Frequency__C?_CP_STA_PPP'',''mat'');',3:4);
c_eval('specEell?.p=specEell?.p(:,[2:end]);',3:4) % take away time row
%c_eval('specEell?.p(specEell?.p<-1) = NaN;',3:4)


% density
c_eval('P?=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'',''mat'');',3:4);
c_eval('scpNe?=c_efw_scp2ne(P?);',3:4);
%c_eval('scpNe?=irf_resamp(scpNe?,diE?);',3:4);
c_eval('peaNe?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''mat'');',3:4);
%c_eval('peaNe?hr=irf_resamp(peaNe?,diE?);',3:4);

% plasma frequency
c_eval('peafpe? = irf_plasma_calc(diB?,peaNe?,peaNe?,peaNe?,peaNe?,''Fpe'');',3:4);
c_eval('scpfpe? = irf_plasma_calc(diB?,irf_tappl(scpNe?,''/20''),scpNe?,scpNe?,scpNe?,''Fpe'');',3:4);



%% make figure
h = irf_plot(6); isub = 1;

if 1 
    hca = h(isub); isub = isub + 1;        
    irf_plot(hca,diE3);
end

if 1 
    hca = h(isub); isub = isub + 1;        
    irf_plot(hca,{peaNe3,irf_tappl(scpNe3,'/20')},'comp');
    irf_legend(hca,{'n_{e,PEA}','n_{e,EFW}/20'},[0.02 0.92])
end

if 1 
    hca = h(isub); isub = isub + 1;        
    irf_spectrogram(hca,specE3); hold(hca,'on');
    set(hca,'yscale','log')
    irf_plot(hca,irf_tappl(peafpe3,'*1e-3'));
    irf_plot(hca,irf_tappl(scpfpe3,'*1e-3'));
    
end

if 1 
    hca = h(isub); isub = isub + 1;        
    irf_spectrogram(hca,specE4); hold(hca,'on');
    set(hca,'yscale','log')
    irf_plot(hca,irf_tappl(peafpe4,'*1e-3'));
    irf_plot(hca,irf_tappl(scpfpe4,'*1e-3'));
end

if 1 
    hca = h(isub); isub = isub + 1;        
    irf_spectrogram(hca,specEpol3,'lin'); hold(hca,'on');   
    set(hca,'yscale','log')
    irf_colormap(hca,'cmap')
    chca = colorbar('peer',hca);
    ylabel(chca,'polarization')
    caxis(hca,[0 1])
end

if 1 
    hca = h(isub); isub = isub + 1;        
    irf_spectrogram(hca,specEell4,'lin'); hold(hca,'on');
    set(hca,'yscale','log')
    irf_colormap(hca,'cmap')    
    chca = colorbar('peer',hca);
    ylabel(chca,'ellipticity')
    caxis(hca,[-1 1])
end

irf_zoom(h,'x',tint)
irf_plot_axis_align


