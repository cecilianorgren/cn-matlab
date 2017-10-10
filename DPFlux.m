figure('name','Event overview')
n_panels=4;
t1=[2007 08 31 10 12 0 0];t2=[2007 08 31 10 26 0 0];
tint2=[toepoch([2007 08 31 10 13 00]) toepoch([2007 08 31 10 24 00])];      
tint=[gsmB3(1,1) gsmB3(end,1)];
h=irf_plot(n_panels);
isub=1;

if 1 % B C3 GSM FGM (1 panel)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,gsmB3fgm);
    ylabel(hca,'B [nT] GSM C3');
    irf_legend(hca,{'x','y','z','|B|'},[0.02 0.05]);
    grid(hca,'off')
end
if 1 % C3_CP_PEA_PITCH_SPIN_DPFlux (1 panels)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DPFlux','sum_dim1','colorbarlabel',{'log_{10} dPF','1/cm^2 s sr keV'},'fitcolorbarlabel');
    set(hca,'yscale','log','tickdir','out');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_e [eV]');
end
if 1 % C3_CP_PEA_3DXPH_DPFlux (1 panels)
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Data__C?_CP_PEA_3DXPH_DPFlux',3);
    [var,dobj,varmat,varunits]=c_caa_var_get(varname); % time, azimuth, polar, energy
    DEF=squeeze(nanmean(nanmean(varmat.data,2),3));
    energy=varmat.dep_x{3}.data(1,:);%+0.5*(varmat.dep_x{3}.df.plus-varmat.dep_x{3}.df.minus);
    en_ind=isnan(energy);
    energy(en_ind)=[];
    DEF(:,en_ind)=[];
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
if 1 % C3_CP_PEA_3DXPH_DPFlux (1 panels)
    hca=h(isub); isub=isub+1;
    varname=irf_ssub('Data__C?_CP_PEA_3DXPH_DPFlux',3);
    res=c_caa_construct_subspin_res_data(varname);
    DEF=squeeze(nanmean(res.data,2)); % sum over pitch angle
    energy=res.en;
    en_ind=isnan(energy);energy(en_ind)=[];
    DEF(:,en_ind)=[];
    
    ud = get(gcf,'userdata');
    t_st_e = double(ud.t_start_epoch);
    pcolor(hca,res.tt-t_st_e,energy,double(log10(DEF)')); 
    
    set(hca,'yscale','log','tickdir','out');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_e [eV]');
    
    hcb = colorbar('peer',hca);    
    irf_colormap(hcb,'space')
    ylabel(hcb,{'log_{10} dPF';varunits},'fontsize',14)
    shading(hca,'flat')
end

set(h(2),'ylim',get(h(3),'ylim'))
for k=2:4; caxis(h(k),[3.5 7.1]); end
irf_colormap('space');
title(h(1),'Event overview: C3');
irf_zoom(h,'x',tint2);

irf_plot_axis_align;
irf_timeaxis(hca,'usefig');
figpos=get(gcf,'position');
set(gcf,'position',[figpos(1) figpos(2) figpos(3)*1.2 figpos(4)]);
set(gcf,'PaperPositionMode','auto'); 