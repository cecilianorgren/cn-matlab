%% Get B,E,P   
c_eval('[caaB?,~,gseB?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');',3:4);
c_eval('[caaE?,~,diE?]=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'');',3:4);
c_eval('[caaP?,~,P?]=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'');',3:4);
c_eval('absB?=irf_abs(gseB?);',3:4);
[caaNi3,~,Ni3]=c_caa_var_get('density__C3_CP_CIS_HIA_ONBOARD_MOMENTS');

%% 1st figure
event=1;
zoom=8;
switch event
    case 1
        date='20070831';
        switch zoom % zoom 20070831 
            case 1
               time='all';
               flow=10;
               tint_zoom=[toepoch([2007 08 31 09 45 00]) toepoch([2007 08 31 12 15 00])]; 
            case 2 % all oscillating area around density jump
               time='h10m15s00-m19s30';
               flow=10;
               tint_zoom=[toepoch([2007 08 31 10 15 00]) toepoch([2007 08 31 10 19 30])];
            case 3 % 
               time='h10m18s30-m19s20';
               flow=10;
               tint_zoom=[toepoch([2007 08 31 10 18 30]) toepoch([2007 08 31 10 19 20])];
            case 4 % two high amp area
               time='h10m18s39-m19s12';
               flow=10;
               tint_zoom=[toepoch([2007 08 31 10 18 39]) toepoch([2007 08 31 10 19 12])];
            case 5 % left high amp area
               time='h10m19s00-12';
               flow=10;
               tint_zoom=[toepoch([2007 08 31 10 19 00]) toepoch([2007 08 31 10 19 12])];
            case 6 % right high amp area
               time='h10m18s39-43';
               flow=10;
               tint_zoom=[toepoch([2007 08 31 10 18 39]) toepoch([2007 08 31 10 18 43])];
            case 7 % right high amp area
               time='h10m19s10-12';
               flow=10;
               tint_zoom=[toepoch([2007 08 31 10 19 10]) toepoch([2007 08 31 10 19 12])];
            case 8 % right high amp area
               time='h10m19s02-08';
               flow=6;
               tint_zoom=[toepoch([2007 08 31 10 19 02]) toepoch([2007 08 31 10 19 08])];
            case 9 % right high amp area
               time='h10m19s05-08';
               flow=10;
               tint_zoom=[toepoch([2007 08 31 10 19 05]) toepoch([2007 08 31 10 19 08])];                           
            case 10 % only highest amplitude area (left)
               time='h10m18s42-43';
               flow=10;
               tint_zoom=[toepoch([2007 08 31 10 18 42]) toepoch([2007 08 31 10 18 43])];
            end
    case 2
        date='20070902';
        switch zoom % zoom 20070902 
            case 1
               time='all';
               flow=1;
               tint_zoom=tint(10,:); 
            case 2
               time='h14m30-36';
               flow=10;
               tint_zoom=[toepoch([2007 09 02 14 30 30]) toepoch([2007 09 02 14 36 00])];
            case 3
               time='';
               flow=10;
               tint_zoom=[toepoch([2007 09 02 14 32 00]) toepoch([2007 09 02 14 33 00])];
            case 4
               time='';
               flow=10;
               tint_zoom=[toepoch([2007 09 02 14 32 20]) toepoch([2007 09 02 14 32 45])];
            case 5 % slightly larger
               time='h14m32s26-32';
               flow=30;
               tint_zoom=[toepoch([2007 09 02 14 32 26]) toepoch([2007 09 02 14 32 32])];
            case 6 % no activity on C3
               time='h14m32s28-29';
               flow=30;
               tint_zoom=[toepoch([2007 09 02 14 32 28]) toepoch([2007 09 02 14 32 29])];   
            case 7 % no activity on C3
               time='h14m32s29-30';
               flow=30;
               tint_zoom=[toepoch([2007 09 02 14 32 29]) toepoch([2007 09 02 14 32 30])];   
case 8 % no activity on C3
               time='hejsan';
               flow=20;
               tint_zoom=[toepoch([2007 09 02 15 46 00]) toepoch([2007 09 02 15 50 00])];   
case 9 % no activity on C3
               time='h14m32s29-30';
               flow=30;
               tint_zoom=[toepoch([2007 09 02 14 32 29]) toepoch([2007 09 02 14 32 30])];   

        end   
    case 3
        date='20070926';
        switch zoom % zoom 20070926 
            case 1
               time='all';
               flow=20;
               tint_zoom=tint(12,:); 
            case 2 % nice waveform
               time='h10m16s34-35';
               flow=15;
               tint_zoom=[toepoch([2007 09 26 10 16 34]) toepoch([2007 09 26 10 16 35])];
            case 3
               time='';
               flow=10;
               tint_zoom=[toepoch([2007 09 26 10 16 00]) toepoch([2007 09 26 10 22 00])];           
            case 5 % some wavy area with calm in the middle
               time='h09m49s05-12';
               flow=10;
               tint_zoom=[toepoch([2007 09 26 09 49 05]) toepoch([2007 09 26 09 49 12])];   
            case 6 % nice waveform Emax~<20
               time='h09m49s07-08';
               flow=20;
               tint_zoom=[toepoch([2007 09 26 09 49 07]) toepoch([2007 09 26 09 49 08])];
            case 7 % no ion data
               time='h10m50s25-26';
               flow=5;
               tint_zoom=[toepoch([2007 09 26 10 50 25]) toepoch([2007 09 26 10 50 26])];
        end
end

%% Calculating flh
tic
i_start=find(diE3(:,1)>tint_zoom(1),1);
i_end=find(diE3(:,1)>tint_zoom(2),1);

i_start_B=find(absB3(:,1)>tint_zoom(1),1,'first');
i_end_B=find(absB3(:,1)<tint_zoom(2),1,'last');

overlap=0;
nfft=1024;
%%
c_eval('diE?zoom=diE?(i_start-nfft:i_end+nfft,:);',3:4);
c_eval('absB?zoom=absB?(i_start_B:i_end_B,:);',3:4);
toc
dt=diE3zoom(2,1)-diE3zoom(1,1);
fs=1/dt;
c_eval('diE?fft=irf_powerfft(diE?zoom,nfft,fs,overlap);',3:4);
%%
%sizeNi3=size(Ni3,1);
%Bindex=fix(linspace(1,size(absB3,1),sizeNi3));
%absB3red=absB3(Bindex,[1 5]);
tic
for t=1:size(absB3zoom,1)
    c_eval('flh?(t,1)=absB?zoom(t,1);',3:4);
    c_eval('flh?(t,2)=irf_plasma_calc(absB?zoom(t,5),0.1,0,2000,3500,''Flh'');',3:4);
end
toc
%% 1st figure
if 1 
    figure
    h=irf_plot(4);
    isub=1;
    diE3zoom=diE3;
    diE4zoom=diE4;
    dt=diE3zoom(2,1)-diE3zoom(1,1);
    fs=1/dt;
    irf_zoom(h,'x',tint_zoom);
    
    if 1 % B GSE C3 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,absB3);
        ylabel(hca,'B [nT]\newline GSE');
        irf_legend(hca,{'x','y','z','|B|','C3'},[0.02 0.05]);
    end
    if 0 % B GSE C4 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,absB4);
        ylabel(hca,'B [nT]\newline GSE');
        irf_legend(hca,{'x','y','z','C4'},[0.02 0.05]);
    end
    if 0 % B ISR2 C3 (1 panel)
        hca=h(isub); isub=isub+1;
        diB3=c_coord_trans('gse','dsi',gseB3,'CL_ID',3);
        irf_plot(hca,diB3);
        ylabel(hca,'B [nT]\newline ISR2');
        irf_legend(hca,{'x','y','z','C3'},[0.02 0.05]);
    end
    if 0 % B ISR2 C4 (1 panel)
        hca=h(isub); isub=isub+1;
        diB4=c_coord_trans('gse','dsi',gseB3,'CL_ID',4);
        irf_plot(hca,diB4);
        ylabel(hca,'B [nT]\newline ISR2');
        irf_legend(hca,{'x','y','z','C4'},[0.02 0.05]);
    end   
    if 0 % |B| (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,absB3(:,[1 5]),'g'); hold(hca,'on');
        irf_plot(hca,absB4(:,[1 5]),'b'); hold(hca,'on');
        ylabel(hca,'|B| [nT]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 1 % E x-comp in ISR2 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,diE3zoom(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,diE4zoom(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X}[mV/m]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 1% E y-comp in ISR2 (1 panel)
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
    if 0 % Ex Ey filtered (2 panels)        
        c_eval('E?_filt=irf_filt(diE?zoom(:,1:3),flow,180,fs,3);',3:4);
        hca=h(isub); isub=isub+1;
        irf_plot(hca,E3_filt(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,E4_filt(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X}[mV/m]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);       
        irf_legend(hca,{['f_{cut}=',num2str(flow),'Hz']},[0.02 0.90]); 
   
        hca=h(isub); isub=isub+1;
        irf_plot(hca,E3_filt(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,E4_filt(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{Y}[mV/m]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);  
        set(hca,'ColorOrder',[0 0 0]);       
        irf_legend(hca,{['f_{cut}=',num2str(flow),'Hz']},[0.02 0.90]); 
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
        irf_spectrogram([hca1 hca2],diE3fft);hold(hca,'on');
        irf_spectrogram([hca3 hca4],diE4fft);hold(hca,'on');
        ylabel(hca1,'f [Hz]');
        ylabel(hca2,'f [Hz]');
        ylabel(hca3,'f [Hz]');
        ylabel(hca4,'f [Hz]');
        irf_legend(hca1,'C3_X',[0.02 0.90]);
        irf_legend(hca2,'C3_Y',[0.02 0.90]);
        irf_legend(hca3,'C4_X',[0.02 0.90]);
        irf_legend(hca4,'C4_Y',[0.02 0.90]);
        
        
        %set(hca1,'yscale','log');
        %set(hca1,'ytick',[1 1e1 1e2 1e3]);
    end
    if 0%i~=1 % E power spectrum with f_lh inlaid (4 panels)
        hca1=h(isub);isub=isub+1;
        hca2=h(isub);isub=isub+1;
        hca3=h(isub);isub=isub+1;
        hca4=h(isub);isub=isub+1; 
        
        %c_eval('diE?fft=irf_powerfft(diE?zoom,nfft,fs,overlap);',3:4);
        irf_spectrogram([hca1 hca2],diE3fft);hold(hca1,'on');hold(hca2,'on');
        irf_spectrogram([hca3 hca4],diE4fft);hold(hca3,'on');hold(hca4,'on');
       
        irf_plot(hca1,flh3);
        irf_plot(hca2,flh3);
        irf_plot(hca3,flh4);
        irf_plot(hca4,flh4);
        
        set(hca1,'yscale','log');
        set(hca1,'ytick',[1 1e1 1e2]);
        set(hca2,'yscale','log');
        set(hca2,'ytick',[1 1e1 1e2]);
        set(hca3,'yscale','log');
        set(hca3,'ytick',[1 1e1 1e2]);
        set(hca4,'yscale','log');
        set(hca4,'ytick',[1 1e1 1e2]);       
        ylabel(hca1,'f [Hz]');
        ylabel(hca2,'f [Hz]');
        ylabel(hca3,'f [Hz]');
        ylabel(hca4,'f [Hz]');
        irf_legend(hca1,'C3_X',[0.02 0.90]);
        irf_legend(hca2,'C3_Y',[0.02 0.90]);
        irf_legend(hca3,'C4_X',[0.02 0.90]);
        irf_legend(hca4,'C4_Y',[0.02 0.90]);
    end
    if 0%DEFlux3==1 % Electron data C3 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([5.9 7.6]);
        set(hca,'yscale','log','ylim',[100 3e4]);
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        irf_legend(hca,'C3',[0.98 0.05],'color','k');
        ylabel(hca,'E_e [eV]');
    end
    if 0%DEFlux4==1 % Electron data C4 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'Data__C4_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([5.9 7.6]);
        set(hca,'yscale','log','ylim',[100 3e4]);
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        irf_legend(hca,'C4',[0.98 0.05],'color','k');
        ylabel(hca,'E_e [eV]');        
    end
    if 0%DPFlux3==1 % Electron data C3 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DPFlux','sum_dim1','colorbarlabel','log_{10} dPF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([5.9 7.6]);
        set(hca,'yscale','log','ylim',[100 3e4]);
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        irf_legend(hca,'C3',[0.98 0.05],'color','k');
        ylabel(hca,'E_e [eV]');
    end
    if 0%DPFlux4==1 % Electron data C4 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'Data__C4_CP_PEA_PITCH_SPIN_DPFlux','sum_dim1','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([5.9 7.6]);
        set(hca,'yscale','log','ylim',[100 3e4]);
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        irf_legend(hca,'C4',[0.98 0.05],'color','k');
        ylabel(hca,'E_e [eV]');        
    end
    if 0 % Ion data C3 (HIA) (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'flux__C3_CP_CIS_HIA_HS_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([3.9 6.1]);
        set(hca,'yscale','log');
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel(hca,'E_i [eV]');
        irf_legend(hca,'C3',[0.98 0.05],'color','k')
    end
    if 0 % Ion data C4 (HIA) (1 panel)
        % is the following right? i changed the variable to plot, but does it
        % still plot an energy?
        hca=h(isub); isub=isub+1; %
        irf_plot(hca,'flux__C4_CP_CIS_HIA_HS_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([3.9 6.1]);
        set(hca,'yscale','log');
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel(hca,'E_i [eV]');
        irf_legend(hca,'C4',[0.98 0.05],'color','k')
        %irf_colormap('default');     
    end  
    if 0 % Ion data C3 (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'flux__C3_CP_CIS_CODIF_H1_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([3.9 6.1]);
        set(hca,'yscale','log');
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel(hca,'E_i [eV]');
        irf_legend(hca,'C3',[0.98 0.05],'color','k')
    end
    if 0 % Ion data C4 (CODIF) (1 panel)
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
    if 0 % Ion density C3 (HIA)
        hca=h(isub); isub=isub+1;
        [caaNi3,~,Ni3]=c_caa_var_get('density__C3_CP_CIS_HIA_ONBOARD_MOMENTS');
        irf_plot(hca,Ni3);
        irf_legend('C3',[0.02 0.93]);
        ylabel(gca,'N_{i} [cm^{-3}]')
    end
    if 0 % Spacecraft potentials (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,{P3,P4},'comp');
        ylabel(hca,'V_{SC} [V]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end   
    if 0 % H+ desnity (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;    
        c_eval('[caaNHi?,~,NHi?]=c_caa_var_get(''density__C3_CP_CIS_CODIF_HS_H1_MOMENTS'');',3:4);    
        irf_plot(hca,NHi3,'g');hold(hca,'on');
        irf_plot(hca,NHi4,'b');hold(hca,'on');
        set(hca,'ColorOrder',[[0 1 0]; [0 0 1]])
        irf_legend({'C3' 'C4'},[0.02 0.93]);
        ylabel(gca,'N_{H^+} [cm^{-3}]')
    end
    if 0 % H+ desnity (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;    
        [caaNHei3,~,NHei3]=c_caa_var_get('density__C3_CP_CIS_CODIF_HS_He1_MOMENTS');
        irf_plot(hca,NHei3);
        ylabel(hca,'N_{He^+} [cm^{-3}]')
    end
    if 0 % O+ desnity (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;    
        [caaNOi3,~,NOi3]=c_caa_var_get('density__C3_CP_CIS_CODIF_HS_O1_MOMENTS');
        irf_plot(hca,NOi3);
        ylabel(hca,'N_{O^+} [cm^{-3}]')
    end
    if 0 % Spacecraft potentials (1 panel)
    hca=h(isub); isub=isub+1;
    c_eval('[caaP?,~,P?]=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'');',3:4);

    irf_plot(hca,P3,'g'); hold(hca,'on');
    irf_plot(hca,P4,'b'); hold(hca,'on');   
   
    ylabel(hca,'V_{SC} [V]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 0 % Ion velocities C3 (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;
        [caaionv3,~,ionv3]=c_caa_var_get('velocity__C3_CP_CIS_CODIF_HS_H1_MOMENTS');
        dionv3=c_coord_trans('gse','dsi',ionv3,'CL_ID',4);     
        irf_plot(hca,dionv3(:,1:3));
        ylabel(hca,'V_i [km/s]')
        irf_legend(hca,'C3',[0.02 0.05])
    end
    if 0 % Ion velocities C4 (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;
        [caaionv4,~,ionv4]=c_caa_var_get('velocity__C4_CP_CIS_CODIF_HS_H1_MOMENTS');
        dionv4=c_coord_trans('gse','dsi',ionv4,'CL_ID',4);
        irf_plot(hca,dionv4(:,1:3));
        ylabel(hca,'V_i [km/s]')
        irf_legend(hca,'C4',[0.02 0.05])
    end
    if 0 % Ion velocities C3 (HIA) (1 panel)
        hca=h(isub); isub=isub+1;
        [caaionv3,~,ionv3]=c_caa_var_get('velocity_isr2__C3_CP_CIS_HIA_ONBOARD_MOMENTS');
        irf_plot(hca,ionv(:,1:3));
        ylabel(hca,'V_i [km/s]')
    end
    if 0 % Ion velocities C4 (HIA) (1 panel)
        hca=h(isub); isub=isub+1;
        [caaionv4,~,ionv4]=c_caa_var_get('velocity_isr2__C4_CP_CIS_HIA_ONBOARD_MOMENTS');
        irf_plot(hca,ionv4(:,1:3));
    end        
    if 0 % ExB velocities C3 (1 panel)
        hca=h(isub); isub=isub+1;
        [caaExB3,~,ExB3]=c_caa_var_get('v_drift_ISR2__C3_CP_EFW_L2_V3D_INERT');
        irf_plot(hca,ExB3(:,1:3));
    end
    if 0 % ExB velocities C4 (1 panel)
        hca=h(isub); isub=isub+1;
        [caaExB4,~,ExB4]=c_caa_var_get('v_drift_ISR2__C4_CP_EFW_L2_V3D_INERT');
        irf_plot(hca,ExB4(:,1:3));
        %irf_plot(hca,'v_drift_ISR2__C3_CP_EFW_L2_V3D_INERT');
    end    
    if 0 % (HIA)-ExB velocities C3 (1 panel)
        hca=h(isub); isub=isub+1;
        newv3=irf_add(1,ionv3,-1,ExB3);
        irf_plot(hca,newv3(:,1:3));
        %irf_plot(hca,'v_drift_ISR2__C3_CP_EFW_L2_V3D_INERT');
    end
    if 0 % (HIA)-ExB velocities C4 (1 panel)
        hca=h(isub); isub=isub+1;
        newv4=irf_add(1,ionv4,-1,ExB4);
        irf_plot(hca,newv4(:,1:3));
    end
    if 0 % ExB xy-comp (2 panels)
        c_eval('[caadiExB?,~,diExB?]=c_caa_var_get(''v_drift_ISR2__C?_CP_EFW_L2_V3D_INERT'');',3:4);
        hca=h(isub); isub=isub+1;        
        irf_plot(hca,diExB3(:,1:2),'g'); hold(hca,'on');
        irf_plot(hca,diExB4(:,1:2),'b'); hold(hca,'on');
        ylabel(hca,'V_{ExB} [km/s]\newline X ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
        
        hca=h(isub); isub=isub+1;        
        irf_plot(hca,diExB3(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,diExB4(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'V_{ExB} [km/s]\newline Y ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end        
    if 0 % Ion velocities xy-comp (CODIF) (2 panels)
        c_eval('[caaionv?,~,ionv?]=c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'');',3:4);
        c_eval('dionv?=c_coord_trans(''gse'',''dsi'',ionv?,''CL_ID'',4);',3:4);   
        
        hca=h(isub); isub=isub+1;
        irf_plot(hca,dionv3(:,1:2),'g'); hold(hca,'on');
        irf_plot(hca,dionv4(:,1:2),'b'); hold(hca,'on');
        ylabel(hca,'V_{i,x} [km/s]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
        
        hca=h(isub); isub=isub+1;
        irf_plot(hca,dionv3(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,dionv4(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'V_{i,y} [km/s]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end    
    if 1 % Electron densities from spacecraft potential (1 panel)
        hca=h(isub); isub=isub+1;
        c_eval('Ne?=c_efw_scp2ne(P?);',3:4);
        
        irf_plot(hca,Ne3,'g'); hold(hca,'on');
        irf_plot(hca,Ne4,'b'); hold(hca,'on');
        ylabel(hca,'N_e [cm^{-3}]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);       
    end
    if 0 % Electron temperature (1 panel)
       hca=h(isub); isub=isub+1; 
       c_eval('[caaTe?par,~,Te?par]=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'');',3:4);
       c_eval('Te?pareV=[Te?par(:,1) Te?par(:,2)*8.61734e-2];',3:4);      
       irf_plot(hca,Te3pareV,'-g'); hold(hca,'on');
       irf_plot(hca,Te4pareV,'b'); hold(hca,'on');
       ylabel(hca,'T_{e,||} [keV]');
       set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
       irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 0 % Ion temperature (HIA) (1 panel)
       hca=h(isub); isub=isub+1; 
       c_eval('[caaTi?,~,Ti?]=c_caa_var_get(''temperature__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',3:4);
       c_eval('Ti?eV=[Ti?(:,1) Ti?(:,2)*8.61734e-2];',3:4);
       ylabel(hca,'T_i [keV]');
       set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
       irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 0 % Ion temperature (CODIF) (1 panel)
       hca=h(isub); isub=isub+1; 
       c_eval('[caaTi?,~,Ti?]=c_caa_var_get(''T__C?_CP_CIS_CODIF_HS_H1_MOMENTS'');',3:4);
       c_eval('Ti?eV=[Ti?(:,1) Ti?(:,2)*8.61734e-2];',3:4);
       ylabel(hca,'T_i [keV]');
       set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
       irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    
    %%
    figure(1)
    irf_zoom(h,'x',tint_zoom)
    irf_plot_axis_align
    add_timeaxis(hca,'usefig');
    position=get(gcf,'Position');
    position(3)=600;
    position(1)=0;
    set(gcf,'PaperPositionMode','auto',...
        'Position',position); 
    eval(['print -dpdf ',date,'_',time,'_exjobb.pdf']);
    
end
%% 2nd figure
if 1 
    figure(2)
    h=irf_plot(12);
    isub=1;
    diE3zoom=diE3;
    diE4zoom=diE4;
    dt=diE3zoom(2,1)-diE3zoom(1,1);
    fs=1/dt;
    irf_zoom(h,'x',tint_zoom);
    
    if 0 % B GSE C3 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,gseB3);
        ylabel(hca,'B [nT]\newline GSE');
        irf_legend(hca,{'x','y','z','C3'},[0.02 0.05]);
    end
    if 0 % B GSE C4 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,gseB4);
        ylabel(hca,'B [nT]\newline GSE');
        irf_legend(hca,{'x','y','z','C4'},[0.02 0.05]);
    end
    if 0 % B ISR2 C3 (1 panel)
        hca=h(isub); isub=isub+1;
        diB3=c_coord_trans('gse','dsi',gseB3,'CL_ID',3);
        irf_plot(hca,diB3);
        ylabel(hca,'B [nT]\newline ISR2');
        irf_legend(hca,{'x','y','z','C3'},[0.02 0.05]);
    end
    if 0 % B ISR2 C4 (1 panel)
        hca=h(isub); isub=isub+1;
        diB4=c_coord_trans('gse','dsi',gseB3,'CL_ID',4);
        irf_plot(hca,diB4);
        ylabel(hca,'B [nT]\newline ISR2');
        irf_legend(hca,{'x','y','z','C4'},[0.02 0.05]);
    end   
    if 1 % |B| (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,absB3(:,[1 5]),'g'); hold(hca,'on');
        irf_plot(hca,absB4(:,[1 5]),'b'); hold(hca,'on');
        ylabel(hca,'|B| [nT]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 1 % E x-comp in ISR2 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,diE3zoom(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,diE4zoom(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X}[mV/m]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 1% E y-comp in ISR2 (1 panel)
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
    if 0 % Ex Ey filtered (2 panels)        
        c_eval('E?_filt=irf_filt(diE?zoom(:,1:3),flow,180,fs,3);',3:4);
        hca=h(isub); isub=isub+1;
        irf_plot(hca,E3_filt(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,E4_filt(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X}[mV/m]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);       
        irf_legend(hca,{['f_{cut}=',num2str(flow),'Hz']},[0.02 0.90]); 
   
        hca=h(isub); isub=isub+1;
        irf_plot(hca,E3_filt(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,E4_filt(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{Y}[mV/m]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);  
        set(hca,'ColorOrder',[0 0 0]);       
        irf_legend(hca,{['f_{cut}=',num2str(flow),'Hz']},[0.02 0.90]); 
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
        irf_spectrogram([hca1 hca2],diE3fft);hold(hca,'on');
        irf_spectrogram([hca3 hca4],diE4fft);hold(hca,'on');
        ylabel(hca1,'f [Hz]');
        ylabel(hca2,'f [Hz]');
        ylabel(hca3,'f [Hz]');
        ylabel(hca4,'f [Hz]');
        irf_legend(hca1,'C3_X',[0.02 0.90]);
        irf_legend(hca2,'C3_Y',[0.02 0.90]);
        irf_legend(hca3,'C4_X',[0.02 0.90]);
        irf_legend(hca4,'C4_Y',[0.02 0.90]);
        
        
        %set(hca1,'yscale','log');
        %set(hca1,'ytick',[1 1e1 1e2 1e3]);
    end
    if 1%i~=1 % E power spectrum with f_lh inlaid (4 panels)
        hca1=h(isub);isub=isub+1;
        hca2=h(isub);isub=isub+1;
        hca3=h(isub);isub=isub+1;
        hca4=h(isub);isub=isub+1; 
        
        %c_eval('diE?fft=irf_powerfft(diE?zoom,nfft,fs,overlap);',3:4);
        irf_spectrogram([hca1 hca2],diE3fft);hold(hca1,'on');hold(hca2,'on');
        irf_spectrogram([hca3 hca4],diE4fft);hold(hca3,'on');hold(hca4,'on');
       
        irf_plot(hca1,flh);
        irf_plot(hca2,flh);
        irf_plot(hca3,flh);
        irf_plot(hca4,flh);
        
        set(hca1,'yscale','log');
        set(hca1,'ytick',[1 1e1 1e2]);
        set(hca2,'yscale','log');
        set(hca2,'ytick',[1 1e1 1e2]);
        set(hca3,'yscale','log');
        set(hca3,'ytick',[1 1e1 1e2]);
        set(hca4,'yscale','log');
        set(hca4,'ytick',[1 1e1 1e2]);       
        ylabel(hca1,'f [Hz]');
        ylabel(hca2,'f [Hz]');
        ylabel(hca3,'f [Hz]');
        ylabel(hca4,'f [Hz]');
        irf_legend(hca1,'C3_X',[0.02 0.90]);
        irf_legend(hca2,'C3_Y',[0.02 0.90]);
        irf_legend(hca3,'C4_X',[0.02 0.90]);
        irf_legend(hca4,'C4_Y',[0.02 0.90]);
    end
    if 0%DEFlux3==1 % Electron data C3 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'Data__C3_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([5.9 7.6]);
        set(hca,'yscale','log','ylim',[100 3e4]);
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        irf_legend(hca,'C3',[0.98 0.05],'color','k');
        ylabel(hca,'E_e [eV]');
    end
    if 0%DEFlux4==1 % Electron data C4 (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'Data__C4_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([5.9 7.6]);
        set(hca,'yscale','log','ylim',[100 3e4]);
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        irf_legend(hca,'C4',[0.98 0.05],'color','k');
        ylabel(hca,'E_e [eV]');        
    end
    if 0 % Ion data C3 (HIA) (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'flux__C3_CP_CIS_HIA_HS_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([3.9 6.1]);
        set(hca,'yscale','log');
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel(hca,'E_i [eV]');
        irf_legend(hca,'C3',[0.98 0.05],'color','k')
    end
    if 0 % Ion data C4 (HIA) (1 panel)
        % is the following right? i changed the variable to plot, but does it
        % still plot an energy?
        hca=h(isub); isub=isub+1; %
        irf_plot(hca,'flux__C4_CP_CIS_HIA_HS_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([3.9 6.1]);
        set(hca,'yscale','log');
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel(hca,'E_i [eV]');
        irf_legend(hca,'C4',[0.98 0.05],'color','k')
        %irf_colormap('default');     
    end  
    if 0 % Ion data C3 (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,'flux__C3_CP_CIS_CODIF_H1_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
        caxis([3.9 6.1]);
        set(hca,'yscale','log');
        set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
        ylabel(hca,'E_i [eV]');
        irf_legend(hca,'C3',[0.98 0.05],'color','k')
    end
    if 0 % Ion data C4 (CODIF) (1 panel)
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
    if 0 % Spacecraft potentials (1 panel)
        hca=h(isub); isub=isub+1;
        irf_plot(hca,{P3,P4},'comp');
        ylabel(hca,'V_{SC} [V]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end   
    if 0 % H+ desnity (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;    
        c_eval('[caaNHi?,~,NHi?]=c_caa_var_get(''density__C3_CP_CIS_CODIF_HS_H1_MOMENTS'');',3:4);    
        irf_plot(hca,NHi3,'g');hold(hca,'on');
        irf_plot(hca,NHi4,'b');hold(hca,'on');
        set(hca,'ColorOrder',[[0 1 0]; [0 0 1]])
        irf_legend({'C3' 'C4'},[0.02 0.93]);
        ylabel(gca,'N_{H^+} [cm^{-3}]')
    end
    if 0 % He+ desnity (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;    
        [caaNHei3,~,NHei3]=c_caa_var_get('density__C3_CP_CIS_CODIF_HS_He1_MOMENTS');
        irf_plot(hca,NHei3);
        ylabel(hca,'N_{He^+} [cm^{-3}]')
    end
    if 0 % O+ desnity (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;    
        [caaNOi3,~,NOi3]=c_caa_var_get('density__C3_CP_CIS_CODIF_HS_O1_MOMENTS');
        irf_plot(hca,NOi3);
        ylabel(hca,'N_{O^+} [cm^{-3}]')
    end
    if 0 % Spacecraft potentials (1 panel)
    hca=h(isub); isub=isub+1;
    c_eval('[caaP?,~,P?]=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'');',3:4);

    irf_plot(hca,P3,'g'); hold(hca,'on');
    irf_plot(hca,P4,'b'); hold(hca,'on');   
   
    ylabel(hca,'V_{SC} [V]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 0 % Ion velocities C3 (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;
        [caaionv3,~,ionv3]=c_caa_var_get('velocity__C3_CP_CIS_CODIF_HS_H1_MOMENTS');
        dionv3=c_coord_trans('gse','dsi',ionv3,'CL_ID',4);     
        irf_plot(hca,dionv3(:,1:3));
        ylabel(hca,'V_i [km/s] \newline C3 CODIF')
        irf_legend(hca,'C3',[0.02 0.05])
    end
    if 1 % Ion velocities C4 (CODIF) (1 panel)
        hca=h(isub); isub=isub+1;
        [caaionv4,~,ionv4]=c_caa_var_get('velocity__C4_CP_CIS_CODIF_HS_H1_MOMENTS');
        dionv4=c_coord_trans('gse','dsi',ionv4,'CL_ID',4);
        irf_plot(hca,dionv4(:,1:3));
        ylabel(hca,'V_i [km/s] \newline C4 CODIF')
        irf_legend(hca,'C4',[0.02 0.05])
    end
    if 1 % Ion velocities C3 (HIA) (1 panel)
        hca=h(isub); isub=isub+1;
        [caaionv3,~,ionv3]=c_caa_var_get('velocity_isr2__C3_CP_CIS_HIA_ONBOARD_MOMENTS');
        irf_plot(hca,ionv3(:,1:3));
        ylabel(hca,'V_i [km/s] \newline C3 HIA')
    end
    if 0 % Ion velocities C4 (HIA) (1 panel)
        hca=h(isub); isub=isub+1;
        [caaionv4,~,ionv4]=c_caa_var_get('velocity_isr2__C4_CP_CIS_HIA_ONBOARD_MOMENTS');
        irf_plot(hca,ionv4(:,1:3));
        ylabel(hca,'V_i [km/s] \newline C4 HIA')
    end        
    if 0 % ExB velocities C3 (1 panel)
        hca=h(isub); isub=isub+1;
        [caaExB3,~,ExB3]=c_caa_var_get('v_drift_ISR2__C3_CP_EFW_L2_V3D_INERT');
        irf_plot(hca,ExB3(:,1:3));
    end
    if 0 % ExB velocities C4 (1 panel)
        hca=h(isub); isub=isub+1;
        [caaExB4,~,ExB4]=c_caa_var_get('v_drift_ISR2__C4_CP_EFW_L2_V3D_INERT');
        irf_plot(hca,ExB4(:,1:3));
        %irf_plot(hca,'v_drift_ISR2__C3_CP_EFW_L2_V3D_INERT');
    end    
    if 0 % (HIA)-ExB velocities C3 (1 panel)
        hca=h(isub); isub=isub+1;
        newv3=irf_add(1,ionv3,-1,ExB3);
        irf_plot(hca,newv3(:,1:3));
        %irf_plot(hca,'v_drift_ISR2__C3_CP_EFW_L2_V3D_INERT');
    end
    if 0 % (HIA)-ExB velocities C4 (1 panel)
        hca=h(isub); isub=isub+1;
        newv4=irf_add(1,ionv4,-1,ExB4);
        irf_plot(hca,newv4(:,1:3));
    end
    if 1 % ExB xy-comp (2 panels)
        c_eval('[caadiExB?,~,diExB?]=c_caa_var_get(''v_drift_ISR2__C?_CP_EFW_L2_V3D_INERT'');',3:4);
        hca=h(isub); isub=isub+1;        
        irf_plot(hca,diExB3(:,1:2),'g'); hold(hca,'on');
        irf_plot(hca,diExB4(:,1:2),'b'); hold(hca,'on');
        ylabel(hca,'V_{ExB} [km/s]\newline X ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
        
        hca=h(isub); isub=isub+1;        
        irf_plot(hca,diExB3(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,diExB4(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'V_{ExB} [km/s]\newline X ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end        
    if 0 % Ion velocities xy-comp (CODIF) (2 panels)
        c_eval('[caaionv?,~,ionv?]=c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'');',3:4);
        c_eval('dionv?=c_coord_trans(''gse'',''dsi'',ionv?,''CL_ID'',4);',3:4);   
        
        hca=h(isub); isub=isub+1;
        irf_plot(hca,dionv3(:,1:2),'g'); hold(hca,'on');
        irf_plot(hca,dionv4(:,1:2),'b'); hold(hca,'on');
        ylabel(hca,'V_{i,x} [km/s]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
        
        hca=h(isub); isub=isub+1;
        irf_plot(hca,dionv3(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,dionv4(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'V_{i,y} [km/s]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end    
    if 0 % Electron densities from spacecraft potential (1 panel)
        hca=h(isub); isub=isub+1;
        c_eval('Ne?=c_efw_scp2ne(P?);',3:4);
        
        irf_plot(hca,Ne3,'g'); hold(hca,'on');
        irf_plot(hca,Ne3,'b'); hold(hca,'on');
        ylabel(hca,'N_e [cm^{-3}]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);       
    end
    if 1 % Electron temperature (1 panel)
       hca=h(isub); isub=isub+1; 
       c_eval('[caaTe?par,~,Te?par]=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'');',3:4);
       c_eval('Te?pareV=[Te?par(:,1) Te?par(:,2)*8.61734e-2];',3:4);      
       irf_plot(hca,Te3pareV,'-g'); hold(hca,'on');
       irf_plot(hca,Te4pareV,'b'); hold(hca,'on');
       ylabel(hca,'T_{e,||} [keV]');
       set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
       irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    if 0 % Ion temperature (1 panel)
       hca=h(isub); isub=isub+1; 
       c_eval('[caaTi?,~,Ti?]=c_caa_var_get(''temperature__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',3:4);
       c_eval('Ti?eV=[Ti?(:,1) Ti?(:,2)*8.61734e-2];',3:4);
       ylabel(hca,'T_i [keV]');
       set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
       irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    end
    
    %%
    irf_zoom(h,'x',tint_zoom)
    irf_plot_axis_align
    add_timeaxis(hca,'usefig');
    position=get(gcf,'Position');
    position(3)=600;
    position(1)=0;
    set(gcf,'PaperPositionMode','auto',...
        'Position',position); 
    eval(['print -dpdf ',date,'_',time,'_4.pdf']);
end