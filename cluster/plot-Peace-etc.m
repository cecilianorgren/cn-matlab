 
%% Figure overview
 cd('/Users/andris/data/JetBraking/Alessandro_20030824_1850_1920');
 tint=[toepoch([2003 8 24 18 50 0]) toepoch([2003 8 24 19 20 0])];
 if 1, % initialize figure
     h=irf_plot(7);
     set(gcf,'defaultlinelinewidth',1.5);
     i_subplot=1;
 end
 if 1, % plot panels
     ic=4;
     if 0, % read FGM data form all sc
         c_eval('[caaB?,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');');
         % c_eval('[caaB?,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_5VPS'');');
         c_eval('B?=irf_abs(B?);');
         c_eval('gsmB?=irf_gse2gsm(B?);');
     end
     if 0, % read CIS HIA/CODIF velocity moments from available s/c
         c_eval('[caaVCIS?,~,VCIS?]=c_caa_var_get(''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',[3]);
         c_eval('[caaVCISH?,~,VCISH?]=c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'');',[4]);
         c_eval('gsmVCIS?=irf_gse2gsm(VCIS?);',[3]);
         c_eval('gsmVCISH?=irf_gse2gsm(VCISH?);',[4]);
     end
     if 1,   % SUBPLOT: FGM Bx
         hca=h(i_subplot);i_subplot=i_subplot+1;
         c_eval('irf_plot(hca,gsmB?(:,[1 2]));',ic);
         ylabel(hca,'B_X [nT] GSM');
         irf_zoom([-30 30],'y',hca);
         irf_legend(hca,{['C' num2str(ic)]},[0.95 0.95],'color','k');
     end
     if 1,   % SUBPLOT: FGM By
         hca=h(i_subplot);i_subplot=i_subplot+1;
         c_eval('irf_plot(hca,gsmB?(:,[1 3]));',ic);
         ylabel(hca,'B_Y [nT] GSM');
         irf_zoom([-30 30],'y',hca);
         irf_legend(hca,{['C' num2str(ic)]},[0.95 0.95],'color','k');
     end
     if 1,   % SUBPLOT: FGM Bz
         hca=h(i_subplot);i_subplot=i_subplot+1;
         c_eval('irf_plot(hca,gsmB?(:,[1 4]));',ic);
         ylabel(hca,'B_Z [nT] GSM');
         irf_zoom([-30 30],'y',hca);
         irf_legend(hca,{['C' num2str(ic)]},[0.95 0.95],'color','k');
     end
     if 1,   % SUBPLOT: CIS Vx
         hca=h(i_subplot);i_subplot=i_subplot+1;
         %        c_eval('irf_plot(hca,gsmVCIS?(:,[1 2]));',ic);   % HIA
         c_eval('irf_plot(hca,gsmVCISH?(:,[1 2]));',ic); % CODIF
         ylabel(hca,'V_X [km/s] GSM');
         irf_zoom([-1500 1500],'y',hca);
         irf_legend(hca,{['C' num2str(ic)]},[0.95 0.95],'color','k');
     end
     if 1,   % SUBPLOT: CIS energy spectrogram
         hca=h(i_subplot); i_subplot=i_subplot+1;
         irf_plot(hca,'flux__C4_CP_CIS_CODIF_H1_1D_PEF','colorbarlabel',...
             'log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
         caxis([3.9 6.1]);
         set(hca,'yscale','log');
         set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
         ylabel(hca,'E [eV]');
         irf_legend(hca,{['C' num2str(ic)]},[0.98 0.05],'color','k')
         irf_colormap('default');
     end
     if 1,    % SUBPLOT: PEACE energy spectrogram
         hca=h(i_subplot); i_subplot=i_subplot+1;
         irf_plot(hca,'Data__C4_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel',...
             'log10 dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
         caxis([5.9 7.6]);
         set(hca,'yscale','log','ylim',[100 3e4]);
         set(hca,'ytick',[1 1e1 1e2 1e3 1e4 1e5])
         irf_legend(hca,{['C' num2str(ic)]},[0.98 0.05],'color','k');
         ylabel('E [eV]');
     end
     if 1,    % SUBPLOT: RAPID electron spectrogram
         hca=h(i_subplot); i_subplot=i_subplot+1;
         irf_plot(hca,'Electron_Dif_flux__C4_CP_RAP_ESPCT6','colorbarlabel',...
             'log10 dF\newline 1/cm^2 s sr keV','fitcolorbarlabel');
         caxis([0.51 4.49]);
         irf_legend(hca,{['C' num2str(ic)]},[0.98 0.85],'color','k')
         set(hca,'yscale','log');
         set(hca,'ytick',[1 1e1 2e1 5e1 1e2 2e2 1e3 1e4 1e5])
     end
     irf_plot_axis_align
     irf_zoom(tint,'x',h);
     irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);
     add_timeaxis(h);
 if 1,         % TMARKS, mark time intervals
     tmarks=toepoch([2003 8 24 18 58 0;2003 8 24 19 03 40]);
     irf_pl_mark(h,tmarks,'black','LineWidth',0.5)
     text_tmarks={'A','B'};
     ypos=ylim(h(1));ypos(2)=ypos(2)+.01*(ypos(2)-ypos(1));ypos(1)=[];
     for j=1:length(tmarks)
         irf_legend(h(1),text_tmarks{j},[tmarks(j)
 ypos],'horizontalalignment','center');
     end
 end
 
 end