if 0 % Largest overview   
figure('name','Event overview')
n_panels=5;
tint=[diE3(1,1) diE3(end,1)];
h=irf_plot(n_panels);
isub=1;
if 1 % B C3 GSM (1 panel)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_abs(gseB3fgm));
    ylabel(hca,'B [nT] GSE');
    irf_legend(hca,{'x','y','z','|B|'},[0.02 0.05]);
end
if 1 % DPFlux3==1 % Electron data C3 (1 panels)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DPFlux','sum_dim1','colorbarlabel','log_{10} dPF\newline #/cm^2 s sr keV','fitcolorbarlabel');
    caxis([5.9 8.6]);
    set(hca,'yscale','log','ylim',[100 3e4]);
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    %irf_legend(hca,'C3',[0.98 0.05],'color','k');
    ylabel(hca,'E_e [eV]');
end
if 1 % Ion data (1 panel)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,'flux__C3_CP_CIS_CODIF_H1_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
    caxis([3.9 6.1]);
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_i [eV]');
    %irf_legend(hca,'C3',[0.98 0.05],'color','k')  
end
if 1 % Spacecraft potentials (1 panel)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,P3);
    ylabel(hca,'V_{SC} [V]');
    %set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    %irf_legend(hca,{'C3','C4'},[0.02 0.05]);
end
if 1 % E C3 ISR2 (1 panel)
    diE3s=diE3;
    diE3s(:,3)=diE3s(:,3)-50;
    hca=h(isub); isub=isub+1;
    irf_plot(hca,diE3s(:,1:3));
    ylabel(hca,'E [mV/m] ISR2');
    irf_legend(hca,{'x','y'},[0.02 0.05]);
end

if 0
    hca=h(isub); isub=isub+1;
    irf_plot(hca,gseExB3); hold on;
    irf_plot(hca,hiaVi3); hold on;
    ylabel(hca,'v_{i,CIS}, v_{ExB} [km/s] GSE')
end
title(h(1),'Event overview C3: 2007-09-02');
irf_zoom(h,'x',tint);
irf_plot_axis_align;
irf_timeaxis(hca,'usefig');

abcde={'a)', 'b)', 'c)', 'd)', 'e)','f)'};
colors={[0 0 0],[0 0 0],[1 1 1],[0 0 0],[0 0 0]};
legloc={[0.02,0.76],[0.02,0.08],[0.02,0.08],[0.02,0.96],[0.02,0.96]};
for k=1:n_panels
    irf_legend(h(k),[abcde(k)],legloc{k},'color',colors{k})
end

set(gcf,'PaperPositionMode','auto');
eval(['print -depsc2 20070902_ov_all.eps']);


end







if 1 % Closer overviewPea

figure('name','Event overview')
n_panels=8;
%tint=[toepoch([2007 08 31 10 18 39.50]) toepoch([2007 08 31 10 18 46.50])]
t1=[2007 08 31 10 12 0 0];t2=[2007 08 31 10 26 0 0];
tint2=[toepoch([2007 08 31 10 13 00]) toepoch([2007 08 31 10 24 00])];
%t1=[2007 9 26 10 45 0 0];
%t2=[2007 9 26 10 55 0 0];
%tint2=[toepoch([2007 9 26 10 45 0 ]) toepoch([2007 9 26 10 54 0 ])];
           
tint=[gsmB3(1,1) gsmB3(end,1)];
h=irf_plot(n_panels);
isub=1;
if 1 %x B C3 GSM FGM (1 panel)
    hca=h(isub); isub=isub+1;
    colors=[0 0 0;1 0 0;0.2 1 0.2];
    set(gca,'colorOrder',colors)
    irf_plot(hca,gsmB4fgm(:,[1 2]),'k');hold(hca,'on');
    irf_plot(hca,gsmB4fgm(:,[1 3]),'r');hold(hca,'on');
    irf_plot(hca,gsmB4fgm(:,[1 4]),'color',[0.2 1 0.2]);hold(hca,'on');
    ylabel(hca,'B [nT] GSM C4');
    irf_legend(gca,{'x','y','z'},[0.98 0.95]);
    set(hca,'ylim',[-5 35])
end
if 0 % deltaB C3 GSM (1 panel)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_add(0.5,gsmB3,-0.5,gsmB4));
    ylabel(hca,'\Delta B [nT] GSM');
    irf_legend(hca,{'x','y','z','|B|'},[0.02 0.05]);
end
if 0 % DPFlux3==1 % Electron data C3 (1 panels)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DPFlux','sum_dim1','colorbarlabel',{'log_{10} dPF','1/cm^2 s sr keV'},'fitcolorbarlabel');
    set(hca,'yscale','log','tickdir','out');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_e [eV]');
end
if 0 % egen DEFlux fr?n 3DXPH_DEFlux % Electron data C3 (1 panels)
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Data__C?_CP_PEA_3DXPH_DEFlux',3);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname); % time, azimuth, polar, energy
    %res=c_caa_construct_subspin_res_data(varname);
    %DEF2=squeeze(nansum(res.data,2));
    DEF=squeeze(nanmean(nanmean(varmat.data,2),3)); % sum over azimuth and polar
    energy=varmat.dep_x{3}.data(1,:);%+0.5*(varmat.dep_x{3}.df.plus-varmat.dep_x{3}.df.minus);
    %DEF=DPF.*repmat(energy,size(DPF,1),1);
    ud = get(gcf,'userdata');
    t_st_e = double(ud.t_start_epoch);
    pcolor(hca,varmat.t-t_st_e,energy,double(log10(DEF)')); 
    
    set(hca,'yscale','log','tickdir','out');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_e [eV]');
    
    hcb = colorbar('peer',hca);  
    set(hcb,'xticklabel',[])
    caxis(hca,[4 7.5])
    ylim=get(hca,'ylim');
    set(hca,'ylim',[100 ylim(2)])
    
    irf_colormap(hcb,'space')
    ylabel(hcb,{'log_{10} dEF';varunits},'fontsize',14)
    shading(hca,'flat')
    irf_timeaxis(hca,'nolabels')
end
if 0 % egen DPFlux % Electron data C3 (1 panels)
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Data__C?_CP_PEA_3DXPH_DPFlux',3);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname); % time, azimuth, polar, energy
    res=c_caa_construct_subspin_res_data(varname);
    DEF2=squeeze(nansum(res.data,2));
    DEF=squeeze(nansum(nansum(varmat.data,2),3));
    energy=varmat.dep_x{3}.data(1,:);%+0.5*(varmat.dep_x{3}.df.plus-varmat.dep_x{3}.df.minus);
    ud = get(gcf,'userdata');
    t_st_e = double(ud.t_start_epoch);
    pcolor(hca,varmat.t-t_st_e,energy,double(log10(DEF)')); 
    
    set(hca,'yscale','log','tickdir','out');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_e [eV]');
    
    hcb = colorbar('peer',hca);    
    irf_colormap(hcb,'space')
    ylabel(hcb,{'log_{10} dPF';varunits},'fontsize',14)
    shading(hca,'flat')
end
if 0 % egen DPFlux % Electron data C3 (1 panels)
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Data__C?_CP_PEA_3DXPH_DPFlux',3);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname); % time, azimuth, polar, energy
    res=c_caa_construct_subspin_res_data(varname);
    DEF2=squeeze(nansum(res.data,2));
    DEF=squeeze(nansum(nansum(varmat.data,2),3));
    energy=varmat.dep_x{3}.data(1,:);%+0.5*(varmat.dep_x{3}.df.plus-varmat.dep_x{3}.df.minus);
    ud = get(gcf,'userdata');
    t_st_e = double(ud.t_start_epoch);
    pcolor(hca,varmat.t-t_st_e,energy,double(log10(DEF)')); 
    
    set(hca,'yscale','log','tickdir','out');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_e [eV]');
    
    hcb = colorbar('peer',hca);    
    irf_colormap(hcb,'space')
    ylabel(hcb,{'log_{10} dPF';varunits},'fontsize',14)
    shading(hca,'flat')
end
if 0 % egen DEFlux fr?n PITCH_SPIN_DPFlux % Electron data C3 (1 panels)
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',3);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname);
    DPF=squeeze(nanmean(varmat.data,2)); % sum over pitch angles
    energy=varmat.dep_x{2}.data(1,:)+0.5*(varmat.dep_x{2}.df.plus-varmat.dep_x{2}.df.minus);
    en_ind=isnan(energy);energy(en_ind)=[];
    DPF(:,en_ind)=[];
    DEF=DPF.*repmat(energy/1000,size(DPF,1),1); % /1000 to change to keV
    
    ud = get(gcf,'userdata');
    t_st_e = double(ud.t_start_epoch);
    pcolor(hca,varmat.t-t_st_e,energy,double(log10(DEF)')); 
    
    set(hca,'yscale','log','tickdir','out');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_e [eV]');
    
    hcb = colorbar('peer',hca);    
    irf_colormap(hcb,'space')
    ylabel(hcb,{'log_{10} dPF';'keV/cm^2-s-str-keV'},'fontsize',14)
    shading(hca,'flat')
end
if 0 % egen DPFlux % Electron data C3 fungerar ej (1 panels)
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Data__C?_CP_PEA_3DXPH_DPFlux',3);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname); % time, azimuth, polar, energy
    
    specrec.p=squeeze(nanmean(nanmean(varmat.data,2),3))';
    specrec.f=varmat.dep_x{3}.data(1,:);%+0.5*(varmat.dep_x{3}.df.plus-varmat.dep_x{3}.df.minus);
    %specrec.df=varmat.dep_x{3}.df;
    
    specrec.t=varmat.t;
    specrec.p_label=['log_{10} dPF' , varunits];
    specrec.f_label='Energy [eV]';
    
    irf_spectrogram(hca,specrec);
    
    %set(hca,'yscale','log','tickdir','out');
    %set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    %ylabel(hca,'E_e [eV]');
    
    %hcb = colorbar('peer',hca);
    
    %irf_colormap(hcb,'space')
    %ylabel(hcb,{'log_{10} dPF';varunits},'fontsize',14)
    %shading(hca,'flat')
end
if 0 % Electron data (1 panel)
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Data__C?_CP_PEA_3DXPH_DEFlux',3);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname); % time, azimuth, polar, energy
    
    DEF=squeeze(nansum(nansum(varmat.data,2),3));
    energy=varmat.dep_x{3}.data(1,:);%+0.5*(varmat.dep_x{3}.df.plus-varmat.dep_x{3}.df.minus);
    ud = get(gcf,'userdata');
    t_st_e = double(ud.t_start_epoch);
    pcolor(hca,varmat.t-t_st_e,energy,double(log10(DEF)')); 
    
    set(hca,'yscale','log','tickdir','out');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_e [eV]');
    
    hcb = colorbar('peer',hca);
    
    irf_colormap(hcb,'space')
    ylabel(hcb,{'log_{10} dEF';varunits},'fontsize',14)
    shading(hca,'flat')
    caxis(hca,[6.5 10.3]);
    
    %irf_legend(hca,'C3',[0.98 0.05],'color','k')  
end
if 1 %x Electron data C4 (1 panel)
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Data__C?_CP_PEA_3DXPH_DEFlux',4);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname); % time, azimuth, polar, energy
    
    DEF=squeeze(nansum(nansum(varmat.data,2),3));
    energy=varmat.dep_x{3}.data(1,:);%+0.5*(varmat.dep_x{3}.df.plus-varmat.dep_x{3}.df.minus);
    ud = get(gcf,'userdata');
    t_st_e = double(ud.t_start_epoch);
    pcolor(hca,varmat.t-t_st_e,energy,double(log10(DEF)')); 
    
    set(hca,'yscale','log','tickdir','out');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_e [eV]');
    
    hcb = colorbar('peer',hca);
    
    irf_colormap(hcb,'space')
    ylabel(hcb,{'log_{10} dEF';varunits},'fontsize',14)
    shading(hca,'flat')
    caxis(hca,[6.5 10.3]);
    hcaylim=get(hca,'ylim')
    set(hca,'ylim',[60 hcaylim(2)])
    irf_legend(hca,'C4',[0.98 0.05],'color','k')  
end
if 0 % make own DEF from DPF % Electron data C3 (1 panels)
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',3);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname);
    dElow=varmat.dep_x{2}.df.minus;
    dEup=varmat.dep_x{2}.df.plus;
    if 0
        energies=[];
    end
    
    if 1
        energies0=double(varmat.dep_x{2}.data); % 2166x44 energy
        energies=zeros(size(energies0));
        for k=2:size(energies0,2)-1
            energies(:,k)=(energies0(:,k-1)+energies0(:,k+1))/2;
        end
        energies(:,1)=energies0(:,1);
        energies(:,end)=energies0(:,end);
    end
    phi=varmat.dep_x{1}.data; % 2166x12 pitch angles
    flux=varmat.data; % particle flux
    time=double(varmat.t);
    flux(isnan(flux))=0.000;
    fluxsum0=sum(flux,2);
    fluxsum=zeros(size(fluxsum0,1),size(fluxsum0,3));

    for k=1:44;
        fluxsum(:,k)=fluxsum0(:,1,k);
    end
    fluxsum(:,36:end)=[];
    energies(:,36:end)=[];
    energies(isnan(energies))=0;
    eng=log10(fluxsum)+log10(energies);
    eng2=10.^eng;
    eng2(isinf(eng2))=9;
    ud = get(gcf,'userdata');
    t_st_e = double(ud.t_start_epoch);
    pcolor(hca,time(:,1)-t_st_e,energies(1,:),double(eng')); 
    shading(hca,'flat')
    set(hca,'ytick',[10^3 10^4]);
    hcb = colorbar('peer',hca);
    ylabel(hcb,{'log_{10} dEF';'keV/cm^2s sr keV'},'fontsize',14)
    %irf_colormap(hca,'space');
    
    set(hca,'yscale','log','TickDir','Out');%'ylim',[100 3e4]);
    %set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    %caxis(hca,[3 6.5])
    %irf_legend(hca,'C3',[0.98 0.05],'color','k');
    ylabel(hca,'E_e [eV]');
end
if 0 % make own DEF from DPF % Electron data C3 (1 panels) version 2
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',3);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname);
    if 1
        energy=varmat.dep_x{2}.data(1,:)+0.5*(varmat.dep_x{2}.df.plus-varmat.dep_x{2}.df.minus);
    
    end
 
    time=double(varmat.t);
    en_nan=isnan(energy);
    energy(en_nan)=[];
    flux=varmat.data;
    flux(:,:,en_nan)=[];
    flux=squeeze(nansum(flux,2)); % sum over phi
    en_flux=double(flux.*repmat(energy,size(flux,1),1));
    ud = get(gcf,'userdata');
    t_st_e = double(ud.t_start_epoch);
    pcolor(hca,time(:,1)-t_st_e,energy(1,:),log10(en_flux')); 
    shading(hca,'flat')
    set(hca,'ytick',[10^3 10^4]);
    hcb = colorbar('peer',hca);
    ylabel(hcb,{'log_{10} dEF';varunits},'fontsize',14)
    irf_colormap(hca,'space');
    
    set(hca,'yscale','log','TickDir','Out');%'ylim',[100 3e4]);
    %set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    %caxis(hca,[3 6.5])
    %irf_legend(hca,'C3',[0.98 0.05],'color','k');
    ylabel(hca,'E_e [eV]');
end
if 1 %x Ion data (1 panel)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,'flux__C4_CP_CIS_CODIF_H1_1D_PEF','sum_dim1','colorbarlabel',{'log_{10} dEF','keV/cm^2 s sr keV'},'fitcolorbarlabel');
    caxis(hca,[1.5 7.5]);
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_i [eV]');
    %irf_legend(hca,'C3',[0.98 0.05],'color','k')  
end
if 0 % Spacecraft potentials (1 panel)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,P3);
    ylabel(hca,'V_{SC} [V]');
    %set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    %irf_legend(hca,{'C3','C4'},[0.02 0.05]);
end
if 0 % Peace electron density (1 panel)    
    hca=h(isub); isub=isub+1;
    irf_plot(hca,peaNe3); hold(hca,'on')
    ylabel(hca,'n_e [cm^{-3}]');
    irf_legend(hca,{'PEACE'},[0.02 0.2]);
end
if 0 % Plasma beta (1 panel)    
    hca=h(isub); isub=isub+1;
    irf_plot(hca,betaEH); hold(hca,'on')
    ylabel(hca,'\beta');
    set(hca,'YScale','log',...
        'ytick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3]);
    %irf_legend(hca,{'PEACE'},[0.02 0.2]);
end
if 1 %x Peace electron density and plasma beta (1 panel)    
    hca=h(isub); isub=isub+1;
    irf_plot(hca,peaNe4);
    set(hca,'box','off')
    hca(2) = axes('Position',get(hca(1),'Position'));
    irf_plot(hca(2),betaEH,'r')
    set(hca(2),'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required
    set(hca(2),'YAxisLocation','right','YScale','log',...
        'ytick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3]);
    ylabel(hca(2),'\beta');
    set(hca(2),'Color','none','box','off'); % color of axis
    set(hca(2),'XColor','k','YColor','k'); % color of axis lines and numbers

    h2=hca(2);
    
    irf_timeaxis(hca(2),'nolabels')
    irf_timeaxis(hca(1),'nolabels')
    
    irf_legend(hca(1),'n_e',[0.62 0.15],'color','k')
    irf_legend(hca(2),'\beta',[0.5 0.62],'color','r')
    ylabel(hca(1),'n_e [cm^{-3}]');
    grid(hca(1),'off');grid(hca(2),'off')
    set(hca(2),'ylim',[0.001 max(betaEH(:,2))*1.1])
    linkaxes(hca,'x')
    %set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    %irf_legend(hca,{'PEACE'},[0.02 0.2]);
    
end
if 0 %x Ion velocity (1 panel)
    hca=h(isub); isub=isub+1;
    colors=[0 0 0;1 0 0;0.2 1 0.2];
    set(hca,'colorOrder',colors)
    irf_plot(hca,gsmhiaVi3(:,[1 2]),'k');hold(hca,'on');
    irf_plot(hca,gsmhiaVi3(:,[1 3]),'r');hold(hca,'on');
    irf_plot(hca,gsmhiaVi3(:,[1 4]),'color',[0.2 1 0.2]);hold(hca,'on');
    ylabel(hca,'V_i [km/s] GSM');
    irf_legend(hca,{'x','y','z'},[0.98 0.95]);
end
if 0 % E C3 ISR2 (1 panel)
    diE3s=diE3; 
    diE3s(:,3)=diE3s(:,3)-50;
    hca=h(isub); isub=isub+1;
    irf_plot(hca,diE3s(:,1:3));hold(hca,'on');
    ylabel(hca,'E [mV/m] ISR2');
    irf_legend(hca,{'x','y'},[0.02 0.05]);
end
if 1 % E C3 ISR2 (1 panel)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,diE3(:,1:3));hold(hca,'on');
    ylabel(hca,'E [mV/m] ISR2');
    irf_legend(hca,{'x','y'},[0.02 0.05]);
end
if 1 % E C4 ISR2 (1 panel)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,diE4(:,1:3));hold(hca,'on');
    ylabel(hca,'E [mV/m] ISR2');
    irf_legend(hca,{'x','y'},[0.02 0.05]);
end
if 0 %x one component of gsmE3 C3 (1 panel)
    % Create lower res data
    
    hca=h(isub); isub=isub+1;
    irf_plot(hca,gsmE4(:,[1 3]),'r');hold(hca,'on');
    irf_plot(hca,gsmE4lowres2(:,[1 3]),'k','linewidth',1.5);
    set(hca,'colororder',[[1 0 0]; [0 0 0]])
    irf_legend(hca,{'450 Hz','0.25 Hz'},[0.98 0.95]);
    ylabel(hca,'E_y [mV/m] GSM');
    %irf_legend(hca,{'x','y'},[0.02 0.05]);
end

if 0 % parallel component of E C3 C4 assuming Epar (1 panel)
    c_eval('diB?=c_coord_trans(''gsm'',''dsi'',gsmB?,''cl_id'',?);',3:4)
    c_eval('diB?=irf_resamp(diB?,cn_toepoch(t1,t2,diE?));',3:4)
    c_eval('[diEt? eb_anglong?]=irf_edb(cn_toepoch(t1,t2,diE?),cn_toepoch(t1,t2,diB?),90,''Epar'');',3:4);
    c_eval('Ep?=irf_dot(diEt?,[diB?(:,1) diB?(:,2:4)./repmat(diB?(:,5),1,3)]);',3:4)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,Ep3(:,[1 2]),'k');hold(hca,'on');
    irf_plot(hca,Ep4(:,[1 2]),'r');
    ylabel(hca,'E_{||} [mV/m] GSM');
    irf_legend(hca,{'C3','C4'},[0.98 0.95]);
end
if 0 % |E] C3 (1 panel)
    absE3=irf_abs(gseE3);
    hca=h(isub); isub=isub+1;
    irf_plot(hca,absE3(:,[1 5]));hold(hca,'on');
    ylabel(hca,'|E| [mV/m]');
    %irf_legend(hca,{'x','y'},[0.02 0.05]);
end
if 0 % Ion velocity from CODIF C3+C4 (3 panels)
    c_eval('[caacodifVi?,~,codifVi?]=c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'');',3:4 );
    c_eval('gsmVi?=c_coord_trans(''gse'',''gsm'',codifVi?,''cl_id'',?);',3:4)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,gsmVi3(:,[1 2]),'k');hold(hca,'on');
    irf_plot(hca,gsmVi4(:,[1 2]),'r');
    ylabel(hca,'v_i_X  GSM');
    irf_legend(hca,{'C3','C4'},[0.98 0.95]);
    hca=h(isub); isub=isub+1;
    irf_plot(hca,gsmVi3(:,[1 3]),'k');hold(hca,'on');
    irf_plot(hca,gsmVi4(:,[1 3]),'r');
    ylabel(hca,'v_i_Y  GSM');
    irf_legend(hca,{'C3','C4'},[0.98 0.95]);
    hca=h(isub); isub=isub+1;
    irf_plot(hca,gsmVi3(:,[1 4]),'k');hold(hca,'on');
    irf_plot(hca,gsmVi4(:,[1 4]),'r');
    ylabel(hca,'v_i_Z  GSM');
    irf_legend(hca,{'C3','C4'},[0.98 0.95]);
end
if 0 % Ion velocity from CODIF+HIA C3+C4 (3 panels)
    c_eval('[caacodifVi?,~,codifVi?]=c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'');',3:4 );
    c_eval('gsmVi?=c_coord_trans(''gse'',''gsm'',codifVi?,''cl_id'',?);',3:4)
    c_eval('[caahiafVi?,~,hiaVi?]=c_caa_var_get(''velocity_gse__C3_CP_CIS_HIA_ONBOARD_MOMENTS'');',3 );
    c_eval('gsmhiaVi?=c_coord_trans(''gse'',''gsm'',codifVi?,''cl_id'',?);',3)
    hca=h(isub); isub=isub+1;
    %irf_plot(hca,gsmVi3(:,[1 2]),'k');hold(hca,'on');
    irf_plot(hca,gsmVi4(:,[1 2]),'r');hold(hca,'on');
    irf_plot(hca,gsmhiaVi3(:,[1 2]),'g');
    ylabel(hca,'v_i_X  GSM');
    irf_legend(hca,{'codifC4','hiaC3'},[0.98 0.95]);
    hca=h(isub); isub=isub+1;
    %irf_plot(hca,gsmVi3(:,[1 3]),'k');hold(hca,'on');
    irf_plot(hca,gsmVi4(:,[1 3]),'r');hold(hca,'on');
    irf_plot(hca,gsmhiaVi3(:,[1 3]),'g');
    ylabel(hca,'v_i_Y  GSM');
    irf_legend(hca,{'codifC4','hiaC3'},[0.98 0.95]);
    hca=h(isub); isub=isub+1;
    %irf_plot(hca,gsmVi3(:,[1 4]),'k');hold(hca,'on');
    irf_plot(hca,gsmVi4(:,[1 4]),'r');hold(hca,'on');
    irf_plot(hca,gsmhiaVi3(:,[1 4]),'g');
    ylabel(hca,'v_i_Z  GSM');
    irf_legend(hca,{'codifC4','hiaC3'},[0.98 0.95]);
end
if 0 % Electron parallel velocity from PEACE C3+C4 (2 panels)
    c_eval('[caaPeaceVpar,~,PeaceVpar?]=c_caa_var_get(''Data_Velocity_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'');',3:4);
    c_eval('pVpar?=c_coord_trans(''gse'',''gsm'',PeaceVpar?,''cl_id'',?);',3:4)
    c_eval('pVpar?=irf_dot(PeaceVpar?,[gsmB?(:,1) gsmB?(:,2:4)./repmat(gsmB?(:,5),1,3)]);',3:4)
    %c_eval('[caaPeaceVper,~,PeaceVper?]=c_caa_var_get(''Data_Velocity_ComponentPerpendicularToMagField__C?_CP_PEA_MOM'');',3:4);
    %c_eval('pVper?=c_coord_trans(''gse'',''gsm'',PeaceVper?,''cl_id'',?);',3:4)
    %c_eval('pVper?=irf_abs(pVper?);',3:4)
    %c_eval('pVpar?=irf_abs(PeaceVpar?);',3:4)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,pVpar3(:,[1 2]),'k');hold(hca,'on');
    irf_plot(hca,pVpar4(:,[1 2]),'r');
    ylabel(hca,'v_{e,||}');
    
    %hca=h(isub); isub=isub+1;
    %irf_plot(hca,pVper3(:,[1 5]),'g');hold(hca,'on');
    %irf_plot(hca,pVper4(:,[1 5]),'b');
    %ylabel(hca,'|v_{e,perp}|');
end
if 1 %x E power spectrum (1 component)
    hca=h(isub);isub=isub+1;

    dt=diE3(2,1)-diE3(1,1);
    fs=1/dt;
    overlap=50;
    nfft=512;
    %c_eval('diE?fft=irf_powerfft(cn_toepoch(t1,t2,diE?),nfft,fs,overlap);',4);
    irf_spectrogram(hca,diE4fft);
    ylabel(hca,'f [Hz]');
    %irf_legend(hca,'y ISR2',[0.02 0.04]);
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3]);
    ylim=get(hca,'ylim');
    k2=colorbar('peer',hca,'location','eastoutside');
    ylabel(k2,'log_{10} (mV/m)^2/Hz','fontsize',12)
    hold(k2,'on')
    %ylabel(hca,'E_y [Hz] C3 ISR2')
    %flhline=zeros(size(gsmB3,1),2);
    %flhline(:,1)=gsmB3(:,1);
    %for k = 1:size(gsmB3,1)
    % flhline(k,2)=irf_plasma_calc(gsmB3(k,5),10,0,1000,2000,'Flh');
    %end
    flhline=irf_plasma_calc(gsmB3fgm(:,[1 5]),10,0,1000,2000,'Flh');
    hold(hca,'on');
    irf_plot(hca,flhline,'linewidth',1);
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3]);
    set(hca,'ylim',ylim)
    eleg=irf_legend(hca,'f_{LH}',[0.02 0.45]);
    caxis(h(6),[-6 1])
end
if 1 %x B power spectrum (1 component)
    hca=h(isub);isub=isub+1;

    dt=gsmB3(2,1)-gsmB3(1,1);
    fs=1/dt;
    overlap=20;
    nfft=512;
    c_eval('gsmB?fft=irf_powerfft(gsmB?,nfft,fs,20);',4);
    irf_spectrogram(hca,gsmB4fft);
    ylabel(hca,'f [Hz]');
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3]);
    ylim=get(hca,'ylim');
    k1=colorbar('peer',hca,'location','eastoutside');
    ylabel(k1,'log_{10} (nT)^2/Hz','fontsize',12)
    hold(k1,'on')
    %irf_legend(hca,'y ISR2',[0.02 0.04]);
    %ylabel(hca,'E_y [Hz] C3 ISR2')
    flhline=zeros(size(gsmB3,1),2);
    flhline(:,1)=gsmB3(:,1);
    if 0
    for k = 1:size(gsmB3,1)
    flhline(k,2)=irf_plasma_calc(gsmB4(k,5),10,0,1000,2000,'Flh');
    end
    else
       flhline(:,2)=irf_plasma_calc(gsmB4(:,5),10,0,1000,2000,'Flh');
    end
    hold(hca,'on');
    irf_plot(hca,flhline,'linewidth',1);
    set(hca,'ylim',ylim);
    bleg=irf_legend(hca,'f_{LH}',[0.02 0.45])
    caxis(hca,[-10 -1])
end
if 0 %x B power spectrum (wavelet) (1 component)
    hca=h(isub);isub=isub+1;

    dt=gsmB3(2,1)-gsmB3(1,1);
    fs=1/dt;
    overlap=20;
    nfft=512;
    c_eval('gsmB?fft=irf_wavelet(gsmB?,''nf'',nfft);',3);
    irf_spectrogram(hca,gsmB3fft);
    ylabel(hca,'f [Hz]');
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3]);
    ylim=get(hca,'ylim');
    k1=colorbar('peer',hca,'location','eastoutside');
    ylabel(k1,'log_{10} (nT)^2/Hz','fontsize',10)
    hold(k1,'on')
    %irf_legend(hca,'y ISR2',[0.02 0.04]);
    %ylabel(hca,'E_y [Hz] C3 ISR2')
    flhline=zeros(size(gsmB3,1),2);
    flhline(:,1)=gsmB3(:,1);
    for k = 1:size(gsmB3,1)
    flhline(k,2)=irf_plasma_calc(gsmB3(k,5),10,0,1000,2000,'Flh');
    end
    hold(hca,'on');
    irf_plot(hca,flhline);
    set(hca,'ylim',ylim);
    bleg=irf_legend(hca,'f_{LH}',[0.02 0.45])
    caxis(h(7),[-10 -1])
end
if 0
    hca=h(isub); isub=isub+1;
    irf_plot(hca,gseExB3);
    ylabel(hca,'v_{ExB} [km/s] GSE')
    irf_legend(hca,{'x','y','z'},[0.02 0.05]);
end
if 0
    hca=h(isub); isub=isub+1;
    irf_plot(hca,codifVi3); hold on;
    ylabel(hca,'v_{i,CIS} [km/s] GSE')
    irf_legend(hca,{'x','y','z'},[0.02 0.05]);
end
abcde={'a)', 'b)', 'c)', 'd)', 'e)','f)','g)','h)'};
colors={[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0],[0 0 0]};
legloc={[0.02,0.91],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96]};
for k=1:n_panels
    irf_legend(h(k),abcde(k),legloc{k},'color',colors{k})
    grid(h(k),'off')
end
    
irf_colormap('space');
title(h(1),'Event overview: C3');
irf_zoom([h h2],'x',tint2);
irf_plot_axis_align(h);
irf_timeaxis(hca,'usefig');

if 1
    %irf_plot_ylabels_align(h);
    irf_timeaxis(hcb,'nolabels')
    irf_plot_axis_align([h h2])
end

figpos=get(gcf,'position');
set(gcf,'position',[figpos(1) figpos(2) figpos(3)*1.2 figpos(4)]);
set(gcf,'PaperPositionMode','auto');
%eval(['print -depsc2 ov_eh.eps']);   



end