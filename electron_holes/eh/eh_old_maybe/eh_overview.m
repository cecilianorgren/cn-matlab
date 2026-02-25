%% load magnetic field
c_eval('caa_load(''C?_CP_FGM_FULL'');',1:4);
c_eval('gseB?=getmat(C?_CP_FGM_FULL,''B_vec_xyz_gse__C?_CP_FGM_FULL'');',1:4);
c_eval('gsmB?=c_coord_trans(''gse'',''gsm'',gseB?,''cl_id'',?);',1:4);
%% load ion velocitites hia
sc=[1 3];
c_eval('caa_load(''C?_CP_CIS-HIA_ONBOARD_MOMENTS'');',sc);
c_eval('gsehiaVi?=getmat(C?_CP_CIS_HIA_ONBOARD_MOMENTS,''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',sc);
c_eval('gsmhiaVi?=c_coord_trans(''gse'',''gsm'',gsehiaVi?,''cl_id'',?);',sc);
%% load ion velocitites codif
sc=[3:4];
c_eval('caa_load(''C?_CP_CIS-HIA_ONBOARD_MOMENTS'');',sc);
c_eval('gsehiaVi?=getmat(C?_CP_CIS_HIA_ONBOARD_MOMENTS,''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',sc);
c_eval('gsmhiaVi?=c_coord_trans(''gse'',''gsm'',gsehiaVi?,''cl_id'',?);',sc);
%% Electric field
if ~exist('diB3','var')
    load diEB
end
%c_eval('diB?=c_coord_trans(''gsm'',''dsi'',gsmB?,''cl_id'',?);',3:4)
%c_eval('diB?=irf_resamp(diB?,cn_toepoch(t1,t2,diE?));',3:4)
c_eval('[diEt? eb_anglong?]=irf_edb(cn_toepoch(t1,t2,diE?),cn_toepoch(t1,t2,diB?),90,''Epar'');',3:4);
c_eval('Ep?=irf_dot(diEt?,[diB?(:,1) diB?(:,2:4)./repmat(diB?(:,5),1,3)]);',3:4)
   
%% Load anisotropy data      
varname=irf_ssub('C?_CP_PEA_PITCH_3DXH_PSD',sc);
%varname=irf_ssub('C?_CP_PEA_3DXPH_DEFlux',sc);

xtra_vn=''; % could be _HEEA or _LEEA when using pitch full
[or_var,dobj,varmat,varunits]=c_caa_var_get(['Data',xtra_vn,'__',varname]);
varname2=irf_ssub(['Sweep_Energy',xtra_vn,'__',varname],sc);
%varname3=irf_ssub('Sweep_PitchAngle__C?_CP_PEA_PITCH_3DXH_DEFlux',sc);
energies=c_caa_var_get(varname2,'mat');
energy_units=c_caa_var_get(varname2,'units'); % These give eV, should be keV
%piichangle=c_caa_var_get(varname3,'mat');

% Taking away all fluxes that are below one, these are something else,
% what?
var=or_var;
var.data(or_var.data == or_var.FILLVAL)=NaN;

% taking mean of asimuthal angles
deflux=squeeze(nanmean(var.data(:,:,:,:),2));

% Making the basic flux structure
b_flux.t=energies(:,1);
b_flux.f=energies(1,2:end); nan_en=~isnan(b_flux.f);
b_flux.f=b_flux.f(nan_en);
b_flux.p=deflux(:,:,nan_en);
b_flux.f_units='keV';

clear ind;
ind{1}=1:2;     % 0-30
ind{2}=3:10;    % 30-150
ind{3}=11:12;   % 150-180

for k=1:3
    flux{k}.t=b_flux.t;
    flux{k}.f=b_flux.f;
    flux{k}.p=squeeze(nanmean(b_flux.p(:,ind{k},:),2));
    flux{k}.p_units=varunits;
    flux{k}.f_units='keV'; % energy_units; % gives eV which is wrong
end

% Huishans recommendation
aind=0;
aind=aind+1;
an{aind}.p=log(flux{2}.p./(flux{3}.p));
an{aind}.t=b_flux.t;
an{aind}.f=b_flux.f;
an{aind}.p_units='log(2*per/\newline (par+par))';
an{aind}.f_units=b_flux.f_units;

%%
%figure(1);
%irf_plot({gsmhiaVi1,gsmhiaVi3},'comp')
%% gse
%h=irf_plot(4);
%h=irf_plot({gseB1,gseB2,gseB3,gseB4},'comp');
%irf_legend(h(1),{'C1','C2','C3','C4'},[0.98 0.98],'color','cluster')
%title(h(1),'B GSE')
%ylabel(h(1),'x');ylabel(h(2),'y');ylabel(h(3),'z');
%c_eval('irf_legend(h(?),{''C1'',''C2'',''C3'',''C4''},[0.98 0.98],''color'',''cluster'');',1:3)
%% gsm
%figure(2);
%h=irf_plot({gsmB1,gsmB2,gsmB3,gsmB4},'comp');
%title(h(1),'B GSM')
%ylabel(h(1),'x');ylabel(h(2),'y');ylabel(h(3),'z');
%c_eval('irf_legend(h(?),{''C1'',''C2'',''C3'',''C4''},[0.98 0.98],''color'',''cluster'');',1:3)

%% plot figure
cd /Users/Cecilia/Data/BM/20070831;
loadData = 0;
sc = 3;
t1 = [2007 08 31 10 12 0.0]; t2 = [2007 08 31 10 26 0.0]; tint = toepoch([t1;t2])';
nPanels = 6;
h = irf_plot(nPanels,'newfigure');
isub=1;
if 1 % B C3 GSM (1 panel)       
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_abs(gsmB3));
    ylabel(hca,'B [nT] GSE');
    irf_legend(hca,{'x','y','z','|B|'},[0.02 0.05]);
    irf_zoom(hca,'y')
end
if 0 % DPFlux3==1 % Electron data C3 (1 panels)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DPFlux','sum_dim1','colorbarlabel','log_{10} dPF\newline #/cm^2 s sr keV','fitcolorbarlabel');
    caxis([5.9 8.6]);
    set(hca,'yscale','log','ylim',[100 3e4]);
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    %irf_legend(hca,'C3',[0.98 0.05],'color','k');
    ylabel(hca,'E_e [eV]');
end
if 0 % Ion data (1 panel)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,'flux__C3_CP_CIS_CODIF_H1_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
    caxis([3.9 6.1]);
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
if 0 % E C3 ISR2 (1 panel)
    diE3s=diE3;
    diE3s(:,3)=diE3s(:,3)-50;
    hca=h(isub); isub=isub+1;
    irf_plot(hca,diE3s(:,1:3));
    ylabel(hca,'E [mV/m] ISR2');
    irf_legend(hca,{'x','y'},[0.02 0.05]);
end    
if 1 % parallel component of E C3 C4 assuming Epar (1 panel)
     hca=h(isub); isub=isub+1;
    irf_plot(hca,Ep3(:,[1 2]),'k');hold(hca,'on');
    irf_plot(hca,Ep4(:,[1 2]),'r');
    ylabel(hca,'E_{||} [mV/m] GSM');
    irf_legend(hca,{'C3','C4'},[0.98 0.95]);
end
if 1 % log(per/(par+apar)/2) , recommended by huishan
    
    hca = h(isub); isub = isub + 1;  
    
    k = 1;
    irf_plot(hca,an{k});
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    hcmap(k)=colorbar('peer',hca);
    ylabel(hca,['Energy \newline [',an{k}.f_units,']']) % ,'fontsize',10
    ylabel(hcmap(k),an{k}.p_units) % ,'fontsize',10
    caxis(hca,[-2 2])
end

% Finishing touches to figure
title(h(1),'Event overview C3: 2007-09-02');
irf_zoom(h,'x',tint);
irf_plot_axis_align;
irf_timeaxis(hca,'usefig');

% Labeling
abcde={'a)', 'b)', 'c)', 'd)', 'e)','f)'};
colors={[0 0 0],[0 0 0],[1 1 1],[0 0 0],[0 0 0]};
legloc={[0.02,0.76],[0.02,0.08],[0.02,0.08],[0.02,0.96],[0.02,0.96]};
for k=1:(isub-1)%nPanels
    irf_legend(h(k),[abcde(k)],legloc{k},'color',colors{k})
    grid(h(k),'off')
end

% Fix ylims for some axes
irf_zoom(h(1),'y')

    
    
