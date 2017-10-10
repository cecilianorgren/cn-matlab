% Making figures for Maria Elena Innocenti's paper
cd /Users/Cecilia/Data/Cluster/20070902/
savePath =  '/Users/Cecilia/Research/LH2/20070902/';

% BM interval
tint = toepoch([2007 09 02 14 20 00;2007 09 02 18 50 00])';

% emi.event20070902_loaddata
% emi.event20070902_plotdata

% intersting time intervals
tint = toepoch([2007 09 02 14 32 34;2007 09 02 14 32 38])';
%    2.0070    0.0090    0.0020    0.0150    0.0470    0.0432
%    2.0070    0.0090    0.0020    0.0150    0.0470    0.0499

%%
sc = 3;
csys = 'gsm';
flim = 0.01;
tool.single_event
doPrint = 0;

di_loc = 9e3*sqrt(n_loc)*sqrt(1/1836); % km
pfftE = irf_powerfft(E,size(E,1),450,0.5);
pfftB = irf_powerfft(B,size(B,1),450,0.5);
EBscale = B0.^2*1e-18/(n_loc*1e6*e*4*pi*1e-7);
 
if 1
    tint_out = tint + 1*[-60 60];
    tint_zoom = tint_out;
    h=irf_plot(6,'newfigure');
    hca = irf_panel('B'); c_eval('irf_plot(hca,irf_tlim(irf_abs(gsmB?),tint_zoom+[-4 4]));',sc); hca.YLabel.String = 'B [nT]'; irf_legend(hca,{'x','y','z'},[0.99 0.95])
    hca = irf_panel('Vi'); c_eval('irf_plot(hca,irf_tlim(gsmVi?,tint_zoom+[-4 4]));',sc); hca.YLabel.String = 'V_i [km/s]'; irf_legend(hca,{'x','y','z'},[0.99 0.95])
    %irf_legend(hca,{['X','Y']},[0.02 0.02])
    %hold(hca,'on')
    %co = [0 0 0; 0 0 1; 1 0 0; 
    % h = ax;        
    % set(h,'ColorOrder',co)
    %c_eval('irf_plot(hca,irf_tlim(gsmVi?,tint_zoom));',sc); hca.YLabel.String = 'V_i [km/s]'; irf_legend(hca,{'x','y','z'},[0.99 0.95])
    hca = irf_panel('E'); c_eval('irf_plot(hca,irf_tlim(gsmE?,tint_zoom+[-4 4]));',sc); hca.YLabel.String = 'E [mv/m]'; irf_legend(hca,{'x','y','z'},[0.99 0.95])
    irf_zoom(h(1:3),'x',tint_zoom)
    irf_pl_mark(h(1:3),tint)
    delete(h(4))
    hca = axes('Position',h(5).Position); delete(h(5)); h(5) = hca;
    irf_match_phibe_vis(hca,'best',phi_E(:,[1 1+i_v]),phi_B,v(i_v),direction,corr_dir,flim,flh_loc);
    
    irf_legend(h(5),['E_n = ' num2str(avEn,'%.f') ' mV/m   E_k = ' num2str(avEk,'%.f') ' mV/m   B_0 = ' num2str(B0,'%.f') ' nT   n = ' num2str(n_loc,'%.1f') ' cc' ],[0.02 -.5])
    irf_legend(h(5),['v_{ExB} = [' num2str(irf_norm(vExB),'%.2f') ']\times' num2str(irf_abs(vExB,1),'%.f') 'km/s' ],[0.02 -.75])
    irf_legend(h(5),['v-v_{ExB} = [' num2str(irf_norm(direction*velocity-vExB),'%.2f') ']\times' num2str(irf_abs(direction*velocity-vExB,1),'%.f') ' km/s' ],[0.02 -1.0])
    irf_legend(h(5),['\theta_{v>v_{ExB}} = ' num2str(acosd(direction*irf_norm(vExB)'),'%.f') ' ^{\circ}'],[0.02 -1.25])
    %text(h(5).XLim(1),h(5).YLim(2)*1,'a')
    delete(h(6))
    irf_plot_zoomin_lines_between_panels(hca,h(5))
    title(h(1),['C' num2str(sc) '  ' upper(csys)])
    irf_zoom(h(1:3),'x',tint_zoom)
    
    if 0
    hca = axes('Position',h(8).Position); delete(h(8)); h(8) = hca;
    %
    v_approx = max(Bz(:,2))*B0*1e-18/mu0/e/max(intEdt(:,[1+i_dir]))/n; % km/s
    % E scale: *1e-6*re_loc^2*1e6
    % B scale: *1e-18*(EBscale^2)
    % EB fit scale 
    Eavhighf = mean(pfftB.p{1}(round(numel(pfftB.p{1})*2/3):end));
    Bavhighf = mean(pfftE.p{1}(round(numel(pfftB.p{1})*2/3):end));
    % plot
    semilogy(hca,pfftE.f/flh_loc,pfftE.p{1},pfftB.f/flh_loc,pfftB.p{1},pfftB.f/flh_loc,pfftB.p{1}*Bavhighf/Eavhighf); xlabel(h(8),'f/f_{LH}')
    %semilogy(hca,pfftE.f*re_loc/v_approx,pfftE.p{1},pfftB.f*re_loc/v_approx,pfftB.p{1}); xlabel(h(8),'f\rho_e/v_{ph,approx}')
    set(hca,'YLim',[1e-12 1e1],'ytick',10.^[-20:5:20])
    irf_legend(hca,{'E','B','B_{shift}'},[0.99 0.95])
    ylabel(hca,'Power')
    end
    %%
    if doPrint; cn.print([ 'C' num2str(sc) '_' csys '_inclov_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_best'],'path',savePath); end
end
if 0 % plot spectra
    %%
    figure(34)
    pfftEk = irf_powerfft(Ek(:,[1 1+i_dir]),size(Ek(:,[1 1+i_dir]),1),450,0.5);
    pfftBz = irf_powerfft(Bz,size(Bz,1),450,0.5);
    Ekavhighf = mean(pfftEk.p{1}(round(numel(pfftEk.p{1})*2/3):end));
    Bzavhighf = mean(pfftBz.p{1}(round(numel(pfftBz.p{1})*2/3):end));
    
    pfftPhiE = irf_powerfft(phi_E(:,[1 1+i_v]),size(phi_E(:,[1 1+i_v]),1),450,0.5);
    pfftIntEdt = irf_powerfft(intEdt(:,[1 1+i_dir]),size(intEdt(:,[1 1+i_dir]),1),450,0.5);
    pfftPhiB = irf_powerfft(phi_B,size(phi_B,1),450,0.5);
    pfftV = pfftIntEdt; pfftV.p = {sqrt((pfftPhiB.p{1})./(pfftIntEdt.p{1}))}; % (km/s)
    
    wavelength = pfftV.p{1}'./pfftV.f;
    krho = 2*pi./wavelength*re_loc;
    
    %wavPhiE = irf_wavelet(phi_E(:,[1 1+i_v]),'Fs',450);
    %wavPhiB = irf_wavelet(phi_B,'Fs',450);
    
    pfftE = irf_powerfft(E,size(E,1),450,0.5);
    pfftB = irf_powerfft(B,size(B,1),450,0.5);
    Bavhighf = mean(pfftB.p{1}(round(numel(pfftB.p{1})*2/3):end));
    Eavhighf = mean(pfftE.p{1}(round(numel(pfftB.p{1})*2/3):end));
    
    ncols = 1;
    nrows = 6;
    isub = 1;
    hca = subplot(nrows,ncols,isub); isub = isub + 1;
    semilogy(hca,pfftE.f/flh_loc,pfftE.p{1},pfftE.f/flh_loc,pfftE.p{2},pfftE.f/flh_loc,pfftE.p{3})
    hca.XLabel.String = 'f/f_{LH}';
    hca.YLabel.String = 'Power [mV^2/m^2Hz]';
    irf_legend(hca,{'X','Y','Z'},[0.99 0.95])
    hca.YLim = [min(pfftE.p{1}) max(pfftE.p{1})].*[0.1 20];
    hca.YTick = 10.^[floor(log10(hca.YLim(1))):2:ceil(log10(hca.YLim(2)))];
    
    tint_str= [datestr(irf_time(tint(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint(2),'epoch>datenum'),'hh:mm:ss.fff')];
    hca.Title.String = tint_str;
    
    
    hca = subplot(nrows,ncols,isub); isub = isub + 1;
    semilogy(hca,pfftB.f/flh_loc,pfftB.p{1},pfftB.f/flh_loc,pfftB.p{2},pfftB.f/flh_loc,pfftB.p{3})
    hca.XLabel.String = 'f/f_{LH}';
    hca.YLabel.String = 'Power [nT^2/Hz]';
    irf_legend(hca,{'X','Y','Z'},[0.99 0.95])
    hca.YLim = [min(pfftB.p{1}) max(pfftB.p{1})].*[0.1 20];
    hca.YTick = 10.^[floor(log10(hca.YLim(1))):2:ceil(log10(hca.YLim(2)))];
    
    hca = subplot(nrows,ncols,isub); isub = isub + 1;
    semilogy(hca,pfftBz.f/flh_loc,pfftBz.p{1},pfftEk.f/flh_loc,pfftEk.p{1})
    hca.XLabel.String = 'f/f_{LH}';
    hca.YLabel.String = 'Power';
    irf_legend(hca,{'B_z','E_k'},[0.99 0.95])
    hca.YLim = [min([pfftBz.p{1},pfftEk.p{1}]) max([pfftBz.p{1},pfftEk.p{1}])].*[0.1 20];
    hca.YTick = 10.^[floor(log10(hca.YLim(1))):2:ceil(log10(hca.YLim(2)))];
    
   
    hca = subplot(nrows,ncols,isub); isub = isub + 1;
    semilogy(hca,pfftPhiB.f/flh_loc,pfftPhiB.p{1},pfftPhiE.f/flh_loc,pfftPhiE.p{1})
    hca.XLabel.String = 'f/f_{LH}';
    hca.YLabel.String = 'Power [V^2/Hz]';
    irf_legend(hca,{['\phi_B (n=' num2str(n_loc,'%.2f') 'cm^{-3}, B_0=' num2str(B0,'%.1f') ')'],['\phi_E (v=' num2str(velocity,'%.f') ' km/s)']},[0.99 0.95])
    hca.YLim = [min(pfftPhiB.p{1}) max(pfftPhiB.p{1})].*[0.1 20];
    hca.YTick = 10.^[floor(log10(hca.YLim(1))):2:ceil(log10(hca.YLim(2)))];
    
    hca = subplot(nrows,ncols,isub); isub = isub + 1;
    %plotyy(hca,pfftPhiB.f/flh_loc,pfftPhiB.p{1},pfftIntEdt.f/flh_loc,pfftIntEdt.p{1}*1e-6,'semilogy','semilogy')
    semilogy(hca,pfftPhiB.f/flh_loc,pfftPhiB.p{1},pfftIntEdt.f/flh_loc,pfftIntEdt.p{1}*1e6)
    hca.XLabel.String = 'f/f_{LH}';
    hca.YLabel.String = 'Power [V^2/Hz]';
    irf_legend(hca,{['\phi_B (n=' num2str(n_loc,'%.2f') 'cm^{-3}, B_0=' num2str(B0,'%.1f') 'nT)'],['\int E dt']},[0.99 0.95])
    hca.YLim = [min([pfftPhiB.p{1},pfftIntEdt.p{1}*1e6]) max([pfftPhiB.p{1},pfftIntEdt.p{1}*1e6])].*[0.1 20];
    hca.YTick = 10.^[floor(log10(hca.YLim(1))):2:ceil(log10(hca.YLim(2)))];
    
    hca = subplot(nrows,ncols,isub); isub = isub + 1;
    plot(hca,pfftV.f/flh_loc,(pfftV.p{1})); hca.YLim = [0 10000];
    semilogy(hca,pfftV.f/flh_loc,(pfftV.p{1})); hca.YLim = [0 10000];
    hca.XLabel.String = 'f/f_{LH}';
    hca.YLabel.String = 'v_{ph} [km/s]';    
    hca.YLim = [min(pfftV.p{1}) max(pfftV.p{1})].*[0.1 1];
    %hca.YLim = [5e1 4e3];
    hca.YTick = 10.^[floor(log10(hca.YLim(1))):1:ceil(log10(hca.YLim(2)))];
    
    
    fig=gcf;
    h = fig.Children;
    linkaxes(h,'x')
    for ii = 2:(numel(h))
        h(ii).XLim = [0 7];  
        hold(h(ii),'on')
        plot(h(ii),flim*[1 1],h(ii).YLim,'color',[0 0 0])
        hold(h(ii),'off')
    end 
    
    %if doPrint; cn.print([ 'C' num2str(sc) '_' csys '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_spectra2'],'path',savePath); end    
    if 0
        %%
        figure(35)
        nrows=2;ncols=1;isub=1;
        hca = subplot(nrows,ncols,isub); isub = isub + 1;
        plot(hca,krho,pfftBz.p{1},'+',krho,pfftEk.p{1},'*')
        hca.XLabel.String = 'k/\rho_{e}';
        hca.YLabel.String = 'Power';
        irf_legend(hca,{'B_z','E_k'},[0.99 0.95])
        hca.YLim = [min([pfftBz.p{1},pfftEk.p{1}]) max([pfftBz.p{1},pfftEk.p{1}])].*[0.1 20];
        hca.YTick = 10.^[floor(log10(hca.YLim(1))):2:ceil(log10(hca.YLim(2)))];
        hca.XLim = [0 5];
        hca.YScale = 'log';
        
        hca = subplot(nrows,ncols,isub); isub = isub + 1;
        plot3(hca,krho,pfftBz.f/flh_loc,pfftBz.p{1},'+',krho,pfftEk.f/flh_loc,pfftEk.p{1},'-')
        hca.XLabel.String = 'k/\rho_{e}';
        hca.YLabel.String = 'f/f_{LH}';
        hca.ZLabel.String = 'Power';
        irf_legend(hca,{'B_z','E_k'},[0.99 0.95])
        %hca.ZLim = [min([pfftBz.p{1},pfftEk.p{1}]) max([pfftBz.p{1},pfftEk.p{1}])].*[0.1 20];
        %hca.ZTick = 10.^[floor(log10(hca.ZLim(1))):2:ceil(log10(hca.ZLim(2)))];
        hca.YLim = [0 5];
        hca.XLim = [0 5];
        hca.ZScale = 'log';
    end
end
 %% Visualize
gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,Ek,En,dEn,dEk,mva_l,mva_v,f_highpass,B0);      
gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B,v,n_loc);

%% Print figures
imwrite(gif_stuff_dir.im,gif_stuff_dir.map,[savePath 'C' num2str(sc) '_' csys '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_dir_aaaaaaa.gif'],'DelayTime',0.01,'LoopCount',inf);
imwrite(gif_stuff_v.im,gif_stuff_v.map,[savePath 'C' num2str(sc) '_' csys '_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_v_.gif'],'DelayTime',0.01,'LoopCount',inf);

%% Construct dispersion relation
% dispersion_relation_compare.m
% makes dispersion relation using phi-matching

%tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 13.0])'; % em waves too in high beta
sc = 3;
csys = 'gsm';
flim = 0.1;
%tool.single_event
%dispersion_relation_yoon2008_run
tool.makedisprel
pfftE = irf_powerfft(E,size(E,1),450,0.5);
pfftB = irf_powerfft(B,size(B,1),450,0.5);

if 0 % Plot constructed dispersion relation
    %% Plot 'dispersion relation'
    figure(33)
    nrows = 5;
    tint_str= [datestr(irf_time(tint(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint(2),'epoch>datenum'),'hh:mm:ss.fff')];
    ax=subplot(nrows,1,1);
    plot(f_center,save_vabs,f_center,irf_abs(save_ExB,1),f_center,save_vabs_ionframe,f_center,save_ExBparv);
    ax(1).XLabel.String = '\omega/\omega_{LH}';
    ax(1).YLabel.String = 'v [km/s]';
    title(ax(1),tint_str)
    irf_legend(ax,{'v_{ph}','v_{ExB}','v_{ph}-v_{ExB}','k\cdot v_{ExB}'},[0.02 0.95])
    
    %ax.YLim = [0 1000];

    subplot(nrows,1,2);
    [ax,p1,p2] = plotyy(f_center,save_cv,f_center,save_c,'semilogy','plot');
    ax(1).XLabel.String = '\omega/\omega_{LH}';
    ax(1).YLabel.String = 'C_v';
    ax(1).YTick = 10.^[-50:50];
    
    ax(1).YLim = [min(save_cv) max(save_cv)].*[0.5 1.5];
    ax(2).XAxisLocation = 'top';
    %ax(1).YLim = [5e1 4e3];
    %ax(1).YTick = 10.^[floor(log10(ax(1).YLim(1))):1:ceil(log10(ax(1).YLim(2)))];
    set(ax(1),'box','off')
    ax(2).YLabel.String = 'C_k';
    %ax(2).YLim = [0 1];
    ax(2).YTick = [0:0.1:1];
    
    ax= subplot(nrows,1,3);
    clim = 0.8;
    plot(ax,save_k(save_c>clim)*re,f_center(save_c>clim),'*',save_k*re,f_center)%,...
           %save_k(save_c>clim)*re,save_f_ionframe(save_c>clim),'*',save_k*re,save_f_ionframe);

    legend(ax(1),['C_k>' num2str(clim)],'location','best')
    ax(1).YLabel.String = '\omega/\omega_{LH}';
    ax(1).XLabel.String = 'k\rho_e';
    if ax(1).XLim(2)<2; ax(1).XLim(2)=2; end    
    
    hca = subplot(nrows,1,4);
    [ax,p1,p2] = plotyy(f_center,[save_Ekmax save_Enmax save_Bmax],f_center,[save_EkoverB save_EnoverB]/units.c,'semilogy','plot');
    ax(1).XLabel.String = '\omega/\omega_{LH}';
    ax(1).YLabel.String = '\delta E_{max} [mV/m], \delta B_{max} [nT]';
    ax(2).YLabel.String = '\delta E_{max}/\delta B_{max}/c';
    ax(1).YLim = [min([save_Ekmax;save_Enmax;save_Bmax]) max([save_Ekmax;save_Enmax;save_Bmax])].*[0.9 1.1];
    ax(1).YTick = 10.^[floor(log10(ax(1).YLim(1))):1:ceil(log10(ax(1).YLim(2)))];
    irf_legend(ax(1),{'E_k','E_n','B','E_k','E_n'},[0.03 0.3])
    
    hca = subplot(nrows,1,5);
    plot(hca,save_k(save_c>clim)*re,save_EkoverB(save_c>clim)/units.c,'*',save_k*re,save_EkoverB/units.c,...
             save_k(save_c>clim)*re,save_EnoverB(save_c>clim)/units.c,'*',save_k*re,save_EnoverB/units.c)
    legend(hca,['C_k>' num2str(clim)],'location','best')
    hca.YLabel.String = '\delta E_{max}/\delta B_{max}/c';
    hca.XLabel.String = 'k\rho_e';
    
    %% Plot 'dispersion relation', same as above basically
    tint_str= [datestr(irf_time(tint(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint(2),'epoch>datenum'),'hh:mm:ss.fff')];
    ax=subplot(3,1,1);
    plot(f_center,save_vabs,f_center,irf_abs(save_ExB,1),f_center,save_vabs_ionframe,f_center,save_ExBparv);
    ax(1).XLabel.String = '\omega/\omega_{LH}';
    ax(1).YLabel.String = 'v [km/s]';
    title(ax(1),tint_str)
    irf_legend(ax,{'v_{ph}','v_{ExB}','v_{ph}-v_{ExB}','k\cdot v_{ExB}'},[0.02 0.95])
    %ax(1).YLim =[0 1000]

    ax= subplot(3,1,3);
    clim = 0.7;
    plot(ax,save_k(save_c>clim)*re,f_center(save_c>clim),'*',save_k*re,f_center)%,...
           %save_k(save_c>clim)*re,save_f_ionframe(save_c>clim),'*',save_k*re,save_f_ionframe);

    legend(ax(1),['C_k>' num2str(clim)])
    ax(1).YLabel.String = '\omega/\omega_{LH}';
    ax(1).XLabel.String = 'k\rho_e';
    if ax(1).XLim(2)<2; ax(1).XLim(2)=2; end

    hold(ax,'off');

    subplot(3,1,2);
    [ax,p1,p2] = plotyy(f_center,save_cv,f_center,save_c,'semilogy','plot');
    ax(1).XLabel.String = '\omega/\omega_{LH}';
    ax(1).YLabel.String = 'C_v';
    ax(2).YLabel.String = 'C_k';
    ax(2).YLim = [0 1];

    %% Plot individual correlation fields
    figure(36)
    flim = 0.01
    tool.single_event
    
    nSave = numel(save_phiB);
    if nSave>19; nSave = 19; end
    h=irf_plot(nSave);
    irf_plot(h(1),{phi_E(:,[1 1+i_v]),phi_B}, 'comp')
    legend(h(1),{['f_{highpass}/f_{LH}=' num2str(flim)],['C=' num2str(save_c(ii),'%.2f')]},'location','bestoutside')
    irf_legend(h(1),{'\phi_E','\phi_B'},[0.98 0.95])
    ylabel(h(1),'\phi [V]')
    for ii = 2:nSave;
        dEmax = max(save_dEk{ii}(:,2));
        dBmax = max(save_dBz{ii}(:,2));
        dEB = dEmax*1e-3/dBmax/1e-9;
        irf_plot(h(ii),{save_phiB{ii},save_phiE{ii}},'comp')
        %irf_legend(h(ii),{['f/f_{LH}=' num2str(f_center(ii)) '+[' num2str(fint(1)) ' ' num2str(fint(2)) ']']},[0.98 0.95])
        %irf_legend(h(ii),{['C=' num2str(save_c(ii),'%.2f')]},[0.02 0.95])
        %irf_legend(h(ii),{['\delta E_{max}/\delta B_{max}/c=' num2str(dEB/units.c,'%.2f')]},[0.02 0.95])
        irf_legend(h(ii),{['E/B/c=' num2str(dEB/units.c,'%.2f')]},[0.02 0.95])
        legend(h(ii),{['f/f_{LH}=' num2str(f_center(ii)) '+[' num2str(fint(1)) ' ' num2str(fint(2)) ']'],['C=' num2str(save_c(ii),'%.2f')]},'location','bestoutside')
        ylabel(h(ii),'\phi')
    end
    tint_str = [datestr(irf_time(tint(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint(2),'epoch>datenum'),'hh:mm:ss.fff')];
    info_str = ['\rho_e = ' num2str(re_loc,'%.1f') ' km, \rho_e = ' num2str(re_loc,'%.1f') ', d_i = ' num2str(di_loc,'%.0f') ' km'];
    title(h(1),[tint_str, '   ' info_str])
    irf_zoom(h,'x',tint)
    irf_plot_axis_align
    irf_plot_ylabels_align
    
    %% Similar as above, but plot dEk and Bz instead
    figure(36)
    flim = 0.01
    tool.single_event
    
    nSave = numel(save_phiB);
    if nSave>19; nSave = 19; end
    h=irf_plot(nSave);
    irf_plot(h(1),{save_dBz{ii},save_dEk{ii}}, 'comp')
    legend(h(1),{['f_{highpass}/f_{LH}=' num2str(flim)],['C=' num2str(save_c(ii),'%.2f')]},'location','bestoutside')
    irf_legend(h(1),{'\delta E_k','\delta B_z'},[0.98 0.95])
    ylabel(h(1),'\phi [V]')
    for ii = 2:nSave;
        nAr = size(save_dBz{ii},1);
        plotyy(h(ii),save_dBz{ii}(:,1),save_dBz{ii}(:,2),save_dEk{ii}(:,1),save_dEk{ii}(:,2))
        %irf_legend(h(ii),{['f/f_{LH}=' num2str(f_center(ii)) '+[' num2str(fint(1)) ' ' num2str(fint(2)) ']']},[0.98 0.95])
        %irf_legend(h(ii),{['C=' num2str(save_c(ii),'%.2f')]},[0.02 0.95])
        %legend(h(ii),{['f/f_{LH}=' num2str(f_center(ii)) '+[' num2str(fint(1)) ' ' num2str(fint(2)) ']'],['C=' num2str(save_c(ii),'%.2f')]},'location','bestoutside')
        ylabel(h(ii),'\phi')
    end
    tint_str = [datestr(irf_time(tint(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint(2),'epoch>datenum'),'hh:mm:ss.fff')];
    info_str = ['\rho_e = ' num2str(re_loc,'%.1f') ' km, d_i = '];
    title(h(1),tint_str)
    irf_zoom(h,'x',tint)
    irf_plot_axis_align
    irf_plot_ylabels_align    
    
    %% Plot velocity
    figure(39); 
    tint_str = [datestr(irf_time(tint(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint(2),'epoch>datenum'),'hh:mm:ss.fff')];
    [ax,p1,p2] = plotyy(f_center,save_v,f_center,save_c);
    ylabel(ax(1),'v [km/s]')
    ylabel(ax(2),'C')
    xlabel('f_{center}/f_{LH}')
    legend('x','y','z')
    title(tint_str)
end

%%
omega_lh = flh_loc*2*pi;
clim = 0.7;
hca = subplot(1,4,1);
plot(hca,save_k(save_c>clim)*re,f_center(save_c>clim),'*',save_k*re,f_center,...
     ky/sqrt(2),x_real_store2/omega_lh,ky/sqrt(2),x_imag_store2/omega_lh)
hca.YLabel.String = '\omega/\omega_{LH}';
hca.XLabel.String = 'k\rho_e';

hca = subplot(1,4,2); 
EBscale = B0.^2/(n*1e6*e*4*pi*1e-7);
semilogx(hca,pfftE.p{1},pfftE.f/(omega_lh/2/pi),pfftB.p{1},pfftB.f/(omega_lh/2/pi))
hca.YLabel.String = '\omega/\omega_{LH}';
hca.XLabel.String = 'Power';
hca.XLim = [1e-9 1e2];


hca = subplot(1,4,3); 
semilogx(hca,pfftE.p{1},pfftE.f/(omega_lh/2/pi),pfftE.p{3},pfftE.f/(omega_lh/2/pi),pfftE.p{3},pfftE.f/(omega_lh/2/pi))
semilogx(hca,pfftB.p{1},pfftB.f/(omega_lh/2/pi),pfftB.p{3},pfftB.f/(omega_lh/2/pi),pfftB.p{3},pfftB.f/(omega_lh/2/pi))
semilogx(hca,sqrt(pfftE.p{1}./pfftB.p{1}),pfftB.f/(omega_lh/2/pi),sqrt(pfftE.p{2}./pfftB.p{2}),pfftB.f/(omega_lh/2/pi),sqrt(pfftE.p{3}./pfftB.p{3}),pfftB.f/(omega_lh/2/pi))
hca.YLabel.String = '\omega/\omega_{LH}';
hca.XLabel.String = 'E/B [10^3 km/s]';
hca.XLim = [1e-1 1e3];

hca = subplot(1,4,4); 
semilogy(hca,pfftE.f/(omega_lh/2/pi),pfftE.p{1},pfftB.f/(omega_lh/2/pi),pfftB.p{1})
tisemilogy(hca,pfftE.f/(omega_lh/2/pi),smooth(pfftE.p{1}),pfftB.f/(omega_lh/2/pi),smooth(pfftB.p{1}))
hca.XLabel.String = '\omega/\omega_{LH}';
hca.YLabel.String = 'Power';
hca.YLim = [1e-9 1e2];

%% Check polarization properties
tint_ov = toepoch([2007 09 02 17 53 00; 2007 09 02 17 56 00])';
tint = tint_ov;
%tint_ov = toepoch([2007 09 02 15 47 00; 2007 09 02 15 48 00])';

polarization = irf_ebsp(irf_tlim(gsmE3,tint+30*[-1 1]),irf_tlim(gsmB3,tint+30*[-1 1]),gsmB3,gsmB3,gsmR3,[2 180],'polarization','fac');

frequency = polarization.f;
time = polarization.t;
Bsum = polarization.bb_xxyyzzss(:,:,4);
Esum = polarization.ee_xxyyzzss(:,:,4);
Esum2D = polarization.ee_ss;
ellipticity = polarization.ellipticity;
dop = polarization.dop;
thetak = polarization.k_tp(:,:,1);
planarity = polarization.planarity;
pfluxz = polarization.pf_xyz(:,:,3)./sqrt(polarization.pf_xyz(:,:,1).^2+polarization.pf_xyz(:,:,2).^2+polarization.pf_xyz(:,:,3).^2);

%%
%tint_ov = toepoch([2007 09 02 15 47 00; 2007 09 02 15 48 00])';

ecfreqs = fce3;
tlimit = tint+30*[-1 1];

% Calculate phase speed v_ph = E/B.
vph = sqrt(Esum./Bsum)*1e6;

Bsumthres = 1e-8;
removepts = find(Bsum < Bsumthres);
ellipticity(removepts) = NaN;
thetak(removepts) = NaN;
dop(removepts) = NaN;
planarity(removepts) = NaN;
vph(removepts) = NaN;
pfluxz(removepts) = NaN;


h=irf_plot(6,'newfigure'); 

hca = irf_panel('B');
irf_plot(hca,irf_tlim(gsmB3,tint+30*[-1 1]))
hca.YLabel.String = 'B [nT] GSM';

hca = irf_panel('E');
irf_plot(hca,irf_tlim(gsmE3,tint+30*[-1 1]))
hca.YLabel.String = 'E [mV/m] GSM';

hca=irf_panel('Esum');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=Esum2D;
    specrec.f_label='';
    specrec.p_label={'log_{10}E^{2}','mV^2 m^{-2} Hz^{-1}'};
    irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  irf_legend(hca,'(a)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
irf_plot(hca,flh3,'linewidth',1.5,'color','w')
hold(hca,'off');
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3]);
caxis(hca,[-6 -1])
ylabel(hca,'f (Hz)','fontsize',12);

hca=irf_panel('Bsum');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=Bsum;
    specrec.f_label='';
    specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
    irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  irf_legend(hca,'(b)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
irf_plot(hca,flh3,'linewidth',1.5,'color','w')
hold(hca,'off');
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3]);
caxis(hca,[-8 -2])
ylabel(hca,'f (Hz)','fontsize',12);

hca=irf_panel('vph');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=vph;
    specrec.f_label='';
    specrec.p_label={'log_{10}E/B','m s^{-1}'};
    irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  irf_legend(hca,'(c)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
irf_plot(hca,flh3,'linewidth',1.5,'color','w')
hold(hca,'off');
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3]);
%caxis(hca,[5 8])
ylabel(hca,'f (Hz)','fontsize',12);

hca=irf_panel('thetak');
  specrec=struct('t',time,'p_label',['V^{2}m^{-2}Hz^{-1}']);
    specrec.f=frequency;
    specrec.p=pfluxz;
    specrec.f_label='';
    specrec.p_label={'S_{z}/|S|'};
    irf_spectrogram(hca,specrec,'lin','donotfitcolorbarlabel');
  irf_legend(hca,'(d)',[0.99 0.98],'color','w','fontsize',12)
  hold(hca,'on');
irf_plot(hca,flh3,'linewidth',1.5,'color','w')
hold(hca,'off');
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3]);
caxis(hca,[-1 1])
ylabel(hca,'f (Hz)','fontsize',12);

  
set(h,'xgrid','off','ygrid','off')
set(h(3:6),'Color',0.7*[1 1 1]);

irf_zoom(h,'x',tlimit);
set(h(1:4),'fontsize',12)

irf_plot_axis_align(h(1:4))

%save plot as .png
%set(gcf, 'InvertHardCopy', 'off');
%set(gcf,'paperpositionmode','auto') % to get the same on paper as on screen
%print('-dpng','-painters','-r600','polarizationeg.png');


