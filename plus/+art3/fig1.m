cd /Users/Cecilia/Data/BM/20070831
t1=[2007 08 31 10 13 00];
t2=[2007 08 31 10 25 00];
tint=toepoch([t1;t2])';


%% plot figure
cd /Users/Cecilia/Data/BM/20070831;
loadData = 0;
sc = 4;
t1 = [2007 08 31 10 12 0.0]; t2 = [2007 08 31 10 26 0.0]; tint = toepoch([t1;t2])';
nPanels = 7;
h = irf_plot(nPanels,'newfigure');
isub=1;
leglocind=cell(nPanels,1);
%irf_colormap('poynting')
labelAV=1;

if 1 % B C3 GSM (1 panel)       
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,irf_abs(gsmB?));',sc)
    if labelAV
        ylabel(hca,'B [nT] GSM');
        irf_legend(hca,{'B_{x}','B_{y}','B_{z}','|B|'},[0.98 0.92]);
    else
        ylabel(hca,'B [nT] ');
        irf_legend(hca,{'x_{GSM}','y_{GSM}','z_{GSM}','|B|'},[0.98 0.92]);
    end
    %irf_legend(hca,{'GSM'},[0.07 0.8]);
    irf_zoom(hca,'y')
    legloccind{isub-1} = [0.02 0.8];
end
if 0 % DPFlux3==1 % Electron data C3 (1 panels)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DPFlux',3),'sum_dim1','colorbarlabel','log_{10} dPF\newline #/cm^2 s sr keV');
    cn.colormap(hca,'default')
    caxis(hca,[4.5 7.6]);
    hold(hca,'on');
    irf_plot(hca,irf_tappl(scP3,'*(-1)'))
    set(hca,'yscale','log','ylim',[83 2.6e4]);
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])    
    %irf_legend(hca,'C3',[0.98 0.05],'color','k');
    ylabel(hca,'E_e [eV]');
    legloccind{isub-1} = [0.02 0.8];
end
if 0 % DEFlux3==1 % Electron data C3 (1 panels)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',sc),'sum_dim1','colorbarlabel','log_{10} dPF\newline #/cm^2 s sr keV','fitcolorbarlabel');
    caxis(hca,[4.5 7.6]);
    hold(hca,'on');
    irf_plot(irf_tappl(scP3,'*(-1)'))
    set(hca,'yscale','log','ylim',[50 3e4]);
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])    
    %irf_legend(hca,'C3',[0.98 0.05],'color','k');
    ylabel(hca,'E_e [eV]');
    legloccind{isub-1} = [0.02 0.9];
end

if 0 % Ion data (1 panel)
    hca=h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('flux__C3_CP_CIS_CODIF_H1_1D_PEF',sc),'colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
    caxis([3.9 6.1]);
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    ylabel(hca,'E_i [eV]');
    %irf_legend(hca,'C3',[0.98 0.05],'color','k')  
end
if 0 % Spacecraft potentials (1 panel)
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,P?)',sc);
    ylabel(hca,'V_{SC} [V]');
    %set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    %irf_legend(hca,{'C3','C4'},[0.02 0.05]);
end
if 1 % Peace electron density (1 panel)    
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,peaNe3);',sc); hold(hca,'on')
    ylabel(hca,'n_e [cm^{-3}]');
    %irf_legend(hca,{''},[0.02 0.2]);
    legloccind{isub-1} = [0.02 0.9];
end
if 1 % E C3 ISR2 (1 panel)
    c_eval('diE?s=diE?; diE?s(:,3)=diE?s(:,3)-50;',sc);
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,diE?s(:,1:3));',sc);
    if labelAV
        ylabel(hca,'E [mV/m] ISR2');
        irf_legend(hca,{'E_{x}','E_{y}'},[0.98 0.92]);
    else
        ylabel(hca,'E [mV/m]');
        irf_legend(hca,{'x_{ISR2}','y_{ISR2}'},[0.98 0.92]);
    end
    %irf_legend(hca,{'ISR2'},[0.03 0.9]);
    %irf_legend(h(k),{'C3','C4'},[0.98 0.95]);
    legloccind{isub-1} = [0.02 0.9];
end 
if 1%1 % HIA ion velocity (1 panel)    
    hca=h(isub); isub=isub+1;
    orgsc=sc;
    sc=3;
    c_eval(' irf_plot(hca,gsmVhia?);',sc); hold(hca,'on')
    if labelAV
        ylabel(hca,'v_{i} [km/s] GSM');
        irf_legend(hca,{'v_{x}','v_{y}','v_{z}'},[0.98 0.92]);
    else    
        ylabel(hca,'v_{i} [km/s]');
        irf_legend(hca,{'x_{GSM}','y_{GSM}','z_{GSM}'},[0.98 0.92]);
    end
    %irf_legend(hca,{'GSM'},[0.05 0.9]);
    legloccind{isub-1} = [0.02 0.9];
    sc=orgsc;
end

  
if 1%1 % log(per/(par+apar)/2) , recommended by huishan
    hca = h(isub); isub = isub + 1;      
    k = 2;       
    anprel = an{k}; anprel.f = anprel.f*1000;
    irf_plot(hca,anprel);
    set(hca,'ylim',[60 2.5e4])
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    hcmap(k)=colorbar('peer',hca);
    ylabel(hca,['E_e  [',an{k}.f_units,']']) % ,'fontsize',10
    ylabel(hcmap(k),an{k}.p_units) % ,'fontsize',10
    ylabel(hcmap(k),'2f_{\perp}/(f_{||+}+f_{||-})') % ,'fontsize',10
    caxis(hca,1.8*[-1 1])
    legloccind{isub-1} = [0.02 0.8];
end
if 1%1 % log(par/apar)
    hca = h(isub); isub = isub + 1;  
    k = 1;
    anprel = an{k}; anprel.f = anprel.f*1000;
    irf_plot(hca,anprel);
    set(hca,'ylim',[60 2.5e4])
    set(hca,'yscale','log');
    set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
    hcmap(k)=colorbar('peer',hca);
    ylabel(hca,['E_e  [',an{k}.f_units,']']) % ,'fontsize',10
    ylabel(hcmap(k),an{k}.p_units) % ,'fontsize',10
    ylabel(hcmap(k),'f_{||+}/f_{||-}')
    caxis(hca,1.8*[-1 1])
    legloccind{isub-1} = [0.02 0.8];
end
if 0%1 % Angle between magnetic field and satelltie spin plane (1 panel)    
    hca=h(isub); isub=isub+1;
    c_eval('irf_plot(hca,eb_anglong?);',sc); hold(hca,'on')
    ylabel(hca,'\Theta_{sp,B} [deg]');
    %irf_legend(hca,{'C3','C4'},[0.02 0.2]);
    legloccind{isub-1} = [0.02 0.9];
end
if 0 % Plasma beta (1 panel)    
    hca=h(isub); isub=isub+1;
    irf_plot(hca,betaEH); hold(hca,'on')
    ylabel(hca,'\beta');
    set(hca,'YScale','log',...
        'ytick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3]);
    %irf_legend(hca,{'PEACE'},[0.02 0.2]);
end
if 0 % Pitch angle distributions
    hca=h(isub); isub=isub+1;
    irf_plot(hca,betaEH); hold(hca,'on')
    ylabel(hca,'\beta');
    set(hca,'YScale','log',...
        'ytick',[1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3]);
    %irf_legend(hca,{'PEACE'},[0.02 0.2]);
end
if 0 % parallel component of E C3 C4 assuming Epar (1 panel)
     hca=h(isub); isub=isub+1;
    irf_plot(hca,EparAC3(:,[1 2]),'k');hold(hca,'on');
    irf_plot(hca,EparAC4(:,[1 2]),'r');
    ylabel(hca,'E_{||} [mV/m] GSM');
    irf_legend(hca,{'C3','C4'},[0.98 0.95]);
end

%cn.colormap(h(3:4),'poynting')
%irf_colormap('poynting')

% Finishing touches to figure
%title(h(1),irf_ssub('Event overview C?',sc));
irf_zoom(h,'x',tint);
irf_plot_axis_align;
irf_timeaxis(hca,'usefig');
irf_pl_mark(h,toepoch([2007 08 31 10 17 30;2007 08 31 10 18 30])')

% Labeling
abcde={'a)', 'b)', 'c)', 'd)', 'e)','f)','g)','h)','i)','j)'};
colors={[0 0 0],[0 0 0],[1 1 1],[0 0 0],[0 0 0]};
%legloc={[0.02,0.76],[0.02,0.8],[0.02,0.08],[0.02,0.08],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96]};

for k=1:(isub-1)%nPanels
    irf_legend(h(k),[abcde(k)],legloccind{k},'color',colors{1})
    grid(h(k),'off')
   % irf_legend(h(k),{'C3','C4'},[0.98 0.95]);
end

% Fix ylims for some axes
irf_zoom(h(1),'y')

    
    
