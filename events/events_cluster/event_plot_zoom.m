%% We have three events, choose one
% 1. 20070831
% 2. 20070902
% 3. 20070926
event=1;
%% Reading in time intervals
[t1,t2,tint_comments]=textread('BM.txt','%s%s%[^\n]');
for j=1:size(t1,1),
    tint(j,1)=iso2epoch(t1{j});
    tint(j,2)=iso2epoch(t2{j});
end
%% Going to the right directory etc
    n_panels0=14;
switch event
    case 1 % 20070831
        date= '20070831';
        file_ext='tmp';
        cd /home/cecilia/data/BM/20070831/
        caa_load
        %%
        tint_zoom=[0 0];
        tint_zoom(1,:)=tint(9,:);
        %tint_zoom(2,:)=[toepoch([2007 08 31 10 15 00]) toepoch([2007 08 31 10 25 00])];
        %tint_zoom(3,:)=[toepoch([2007 08 31 10 18 30]) toepoch([2007 08 31 10 19 30])];
        %tint_zoom(4,:)=[toepoch([2007 08 31 10 18 39]) toepoch([2007 08 31 10 18 48])];
        tint_zoom(2,:)=[toepoch([2007 08 31 10 18 41]) toepoch([2007 08 31 10 18 44])];
        %tint_zoom(6,:)=[toepoch([2007 08 31 10 19 00]) toepoch([2007 08 31 10 19 13])];
        
        DEFlux3=0;
        DEFlux4=0;
        n_panels0=n_panels0-2;
    case 2 % 20070902
        date= '20070902';
        file_ext='tmp';
        cd /home/cecilia/data/BM/20070902/
        caa_load
        tint_zoom=[0 0];
        tint_zoom(1,:)=tint(10,:);
        tint_zoom(2,:)=[toepoch([2007 09 02 14 30 30]) toepoch([2007 09 02 14 36 00])];
        tint_zoom(3,:)=[toepoch([2007 09 02 14 32 00]) toepoch([2007 09 02 14 33 00])];
        tint_zoom(4,:)=[toepoch([2007 09 02 14 32 20]) toepoch([2007 09 02 14 32 45])];
        tint_zoom(5,:)=[toepoch([2007 09 02 14 32 26]) toepoch([2007 09 02 14 32 32])];
        DEFlux3=1;
        DEFlux4=1;
    case 3 % 20070926
        date= '20070926';
        file_ext='tmp';
        
        cd /home/cecilia/data/BM/20070926/
        caa_load
        tint_zoom=[0 0];
        tint_zoom(1,:)=tint(12,:);
%        tint_zoom(2,:)=[toepoch([2007 09 26 09 45 00]) toepoch([2007 09 26 10 55 00])];
        tint_zoom(2,:)=[toepoch([2007 09 26 10 16 00]) toepoch([2007 09 26 10 22 00])];
        tint_zoom(3,:)=[toepoch([2007 09 26 09 45 00]) toepoch([2007 09 26 10 55 00])];
        tint_zoom(4,:)=[toepoch([2007 09 26 10 50 25]) toepoch([2007 09 26 10 50 27])];
%         tint_zoom(4,:)=[toepoch([2007 09 02 14 32 20]) toepoch([2007 09 02 14 32 45])];
%         tint_zoom(5,:)=[toepoch([2007 09 02 14 32 26]) toepoch([2007 09
%         02 14 32 32])];        
        DEFlux3=1;
        DEFlux4=0;    
        n_panels0=n_panels0-1;
end
%% Loading into variables
c_eval('[caaB?,~,gseB?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');',3:4);
%c_eval('[caaE?,~,gseE?]=c_caa_var_get(''E_Vec_xyz_GSE__C?_CP_EFW_L2_E3D_GSE'');',3:4);
c_eval('[caaE?,~,diE?]=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'');',3:4);
c_eval('[caaP?,~,P?]=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'');',3:4);
%c_eval('[caaElec?,~,Elec?]=c_caa_var_get(''Data__C3_CP_PEA_PITCH_SPIN_DEFlux'');',3:4);
%% Stuff (not used)
if 0 % E power spectrum
    dt=diE3(2,1)-diE3(1,1);
    fs=1/dt;
    overlap=0;
    c_eval('diE?fft=caa_powerfft(diE?,1024,fs,overlap)',3:4)
end
%% Loop over all relevant time intervals
for i=1:size(tint_zoom,1)-1
    %%
    if i==1 % Only plot spectrogram for zoomed in intervals
        diE3zoom=diE3;
        diE4zoom=diE4;
    else
        i_start=find(diE3(:,1)>tint_zoom(i,1),1);
        i_end=find(diE3(:,1)>tint_zoom(i,2),1);
        diE3zoom=diE3(i_start:i_end,:);
        diE4zoom=diE4(i_start:i_end,:);
    end
    
    % Plot
    eval(['figure(',num2str(i),')']); 
    n_panels=n_panels0;
    if i==1
        n_panels=n_panels-4;
    end

    h=irf_plot(n_panels);
    isub=1;

    if 1 % B C3 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,gseB3);
        ylabel(hca,'B [nT]\newline GSE');
        irf_legend(hca,{'x','y','z','C3'},[0.02 0.05]);
    end
    if 1 % B C4 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,gseB4);
        ylabel(hca,'B [nT]\newline GSE');
        irf_legend(hca,{'x','y','z','C4'},[0.02 0.05]);
    end
    if 1 % |B| (1 panel)
        hca=h(isub); isub=isub+1;
        c_eval('absB?=irf_abs(gseB?);',3:4);
        irf_plot(hca,absB3(:,[1 5]),'g'); hold(hca,'on');
        irf_plot(hca,absB4(:,[1 5]),'b'); hold(hca,'on');
        ylabel(hca,'|B| [nT]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 1 % E x-comp in ISR2 (1 panels)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,diE3zoom(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,diE4zoom(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X}[mV/m]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 1% E y-comp in ISR2 (1 panels)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,diE3zoom(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,diE4zoom(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{Y}[mV/m]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);       
    end
    if 0 % E x+y-comp in GSE (2 panels)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,gseE3(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,gseE4(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{x,ISR2} [mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);

        hca=h(isub); isub=isub+1;
        irf_plot(hca,gseE3(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,gseE4(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{y,GSE} [mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 0%i~=1 % E power spectrum (4 panels)
        hca1=h(isub);isub=isub+1;
        hca2=h(isub);isub=isub+1;
        hca3=h(isub);isub=isub+1;
        hca4=h(isub);isub=isub+1; 

        dt=diE3zoom(2,1)-diE3zoom(1,1);
        fs=1/dt;
        overlap=0;
        nfft=256;
        c_eval('diE?fft=irf_powerfft(diE?zoom,nfft,fs,overlap);',3:4);
        irf_spectrogram([hca1 hca2],diE3fft);
        irf_spectrogram([hca3 hca4],diE4fft); 
        ylabel(hca1,'f [Hz]');
        ylabel(hca2,'f [Hz]');
        ylabel(hca3,'f [Hz]');
        ylabel(hca4,'f [Hz]');
        irf_legend(hca1,'C3_X',[0.02 0.90]);
        irf_legend(hca2,'C3_Y',[0.02 0.90]);
        irf_legend(hca3,'C3_X',[0.02 0.90]);
        irf_legend(hca4,'C3_Y',[0.02 0.90]);
        %set(hca1,'yscale','log');
        %set(hca1,'ytick',[1 1e1 1e2 1e3]);
    end
    if DEFlux3==1 % Electron data C3 (1 panels)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([5.9 7.6]);
        set(hca,'yscale','log','ylim',[100 3e4]);
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        irf_legend(hca,'C3',[0.98 0.05],'color','k');
        ylabel(hca,'E_e [eV]');
    end
    if DEFlux4==1 % Electron data C4 (1 panels)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'Data__C4_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([5.9 7.6]);
        set(hca,'yscale','log','ylim',[100 3e4]);
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        irf_legend(hca,'C4',[0.98 0.05],'color','k');
        ylabel(hca,'E_e [eV]');        
    end
    if 1 % Ion data (2 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'flux__C3_CP_CIS_CODIF_H1_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([3.9 6.1]);
        set(hca,'yscale','log');
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel(hca,'E_i [eV]');
        irf_legend(hca,'C3',[0.98 0.05],'color','k')


        % is the following right? i changed the variable to plot, but does it
        % still plot an energy?
        hca=h(isub); isub=isub+1; %
        irf_plot(hca,'flux__C4_CP_CIS_CODIF_H1_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([3.9 6.1]);
        set(hca,'yscale','log');
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel(hca,'E_i [eV]');
        irf_legend(hca,'C4',[0.98 0.05],'color','k')
        %irf_colormap('default');     
    end  
    if 1 % Spacecraft potentials (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,{P3,P4},'comp');
        ylabel(hca,'V_{SC} [V]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end

    irf_zoom(h,'x',tint_zoom(i,:));
    irf_plot_axis_align
    add_timeaxis(hca,'usefig');
    position=get(gcf,'Position');
    position(3)=500;
    position(1)=0;
    set(gcf,'PaperPositionMode','auto',...
        'Position',position); 
    eval(['print -dpdf ',date,'_',num2str(i),'_',file_ext,'_hej.pdf']);

    if 0
        figure(2)
    end
end
%%
for l=1:fix(size(absB3,1)/1000):size(absB3,1)
    if 1
         F_lh(l) = irf_plasma_calc(absB3(l,2),0.1,0,1000,1000,'Flh');
    end
end