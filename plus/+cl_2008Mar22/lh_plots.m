% Figure for for Maria Elenas paper
cd /Users/Cecilia/Data/Cluster/20070902/
savePath = '/Users/Cecilia/Research/LH2/20070902/';

% Time intervals
tint = toepoch([2007 09 02 15 47 43.2; 2007 09 02 15 47 49.9])';
tint = toepoch([2007 09 02 15 47 45.3; 2007 09 02 15 47 46.7])';
tint_ov = toepoch([2007 09 02 15 46 30; 2007 09 02 15 48 50])';
tint = toepoch([2008 04 22 18 07 33; 2008 04 22 18 07 36.9])';
% Correlation
sc = 3;
csys = 'gsm';

c_eval('R_loc = irf_resamp(gsmR?,tint(1),''nearest'');',sc); R_loc = R_loc(2:4);
c_eval('Vi_loc = irf_resamp(gsmVi?,tint(1),''nearest'');',sc); Vi_loc = Vi_loc(2:4);
c_eval('Ti_loc = irf_resamp(Ti?,tint(1),''nearest'');',sc); Ti_loc = Ti_loc(2);
c_eval('Te_loc = irf_resamp(parTe?,tint(1),''nearest'');',sc); Te_loc = Te_loc(2);


% phi colors
colors_phi = [0.7500    0.1250    0.0980; 
              0.9290    0.6940    0.1250
              0.8500    0.3250    0.0980; 
              0 0 1; 
              1    0.3250    0.0980;   
              0.9290    0.6940    0.1250
              0.8500    0.3250    0.0980; 
              0.9290    0.6940    0.1250];
% length colors              
lengthcolors = [0    0.4470    0.7410; 
                0.4660    0.6740    0.1880; 
                0.4940    0.1840    0.5560;
                0.3010    0.7450    0.9330;
                0.8500    0.3250    0.0980; 
                0.9290    0.6940    0.1250; 
                1 0 0];
                
doPrint = 0;
%% Waveform figure
if 1 % Make main figure
    tint = toepoch([2008 04 22 18 07 33; 2008 04 22 18 07 34])';
    tint_out = tint + 1*[-60 60];
    tint_zoom = tint_out;
    h = irf_plot(6,'newfigure');
    % B
    hca = irf_panel('B'); c_eval('irf_plot(hca,irf_tlim(irf_abs(gsmB?),tint_zoom+[-4 4]));',sc); 
    hca.YLabel.String = 'B [nT]'; irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.99 0.95])
    grid(hca,'off'); 
    %hca.YLim = [0 50];
    hca.Title.String = 'Cluster 3';
    %if hca.YLim(1)>0; hca.YLim = hca.YLim.*[0 1]; end; hca.YTick = [10:10:70];
    irf_legend(hca,'a)',[0.02 0.95]);
    
    % n
    %hca = irf_panel('n'); c_eval('irf_plot(hca,{irf_tlim(peaNe?,tint_zoom+[-8 8]),irf_tlim([scpNe?(:,1) scpNe?(:,2)*0.1],tint_zoom+[-4 4])},''comp'');',sc); 
    %hca.YLabel.String = 'n_e [cm^{-3}]'; irf_legend(hca,{'PEACE','EFW'},[0.99 0.95])
    %grid(hca,'off'); 
    %hca.YTick = [0:0.2:10];
    %hca.YLim = [0 50];
    
    % Vi
    hca = irf_panel('Vi'); c_eval('irf_plot(hca,irf_tlim(gsmVi?,tint_zoom+[-4 4]));',sc); 
    %hca = irf_panel('Vi'); c_eval('irf_plot(hca,{irf_tlim(gsmExB?,tint_zoom+[-4 4]),irf_tlim(gsmVi?,tint_zoom+[-4 4])},''comp'');',sc); 
    hca.YLabel.String = 'V_i [km/s]'; irf_legend(hca,{'V_x','V_y','V_z'},[0.99 0.95])
    grid(hca,'off')
    irf_legend(hca,'b)',[0.02 0.95]);
    % E
    hca = irf_panel('E'); c_eval('irf_plot(hca,irf_tlim(gsmE?,tint_zoom+[-4 4]));',sc); 
    hca.YLabel.String = 'E [mv/m]'; irf_legend(hca,{'E_x','E_y','E_z'},[0.99 0.95])
    grid(hca,'off')
    hca.YLim = [-50 50];
    hca.YTick = [-40:20:40]; 
    irf_legend(hca,'c)',[0.02 0.95]);
    
    irf_zoom(h(1:4),'x',tint_zoom)
    irf_pl_mark(h(1:4),tint)
    hca = irf_panel('delete for space'); %delete(hca)    
    
    grid(hca,'off')    
    irf_timeaxis(h(3))
    %colors_phi = [0 0 0; 1 0 0];
    
    phi_labels_all = {{'\Phi_{\delta E_{k}}','\Phi_{\delta B_{||}}'},...
                      {'\Phi_{\delta E{k}}','\Phi_{\delta B{||}}'},...
                      {'\Phi_{\delta E_{\perp}}','\Phi_{\delta B_{||}}'},...
                      {'\Phi_{\delta E{\perp}}','\Phi_{\delta B{||}}'},...
                      {'\Phi_{\delta E}','\Phi_{\delta B}'}};
    phi_labels = phi_labels_all{5};
                    
    
    % low f phi
    flim = 0.1;
    tool.single_event       
    hca = irf_panel('Phi: low f');
    hca.ColorOrder = colors_phi;
    % move panel slightly upwards to make room for km labels on panel below
    hca.Position = hca.Position + [0 0.02 0 0];
    irf_plot(hca,{phi_E(:,[1 1+i_v]),phi_B},'comp')    
    hca.YLabel.String = '\Phi [V]'; 
    %hleg = irf_legend(hca,{'\Phi_E','\Phi_B'},[0.99 0.95]);
    hleg = irf_legend(hca,phi_labels,[0.99 0.95]);
    %hleg = irf_legend(hca,{'\Phi_{\delta E{\perp}}','\Phi_{\delta B{||}}'},[0.99 0.95]);
    %hleg = irf_legend(hca,{'\Phi_{\delta E}','\Phi_{\delta B}'},[0.99 0.95]);
    hleg(1).Color = colors_phi(1,:);
    hleg(2).Color = colors_phi(2,:);
    %irf_legend(hca,{['v_{ph}=' num2str(velocity,'%.0f') ' km/s']},[0.02 0.95])
    irf_legend(hca,{['\omega > ' num2str(flim,'%.2f') ' \Omega_{LH}']},[0.02 0.1])
    irf_zoom(hca,'x',tint)
    irf_zoom(hca,'y',[-25 25])    
    axtop = add_length_on_top(hca,velocity,1);
    axtop.XLabel.String = 'Length [km]';    
    % add B labels on right
    yratio = phi_B(100,2)/Bz(100,2);            
    axtop.YLim = hca.YLim/yratio;
    axtop.YLabel.String = '\delta B_{||} [nT]'; 
    axtop.YTickLabelMode = 'auto';
    axtop.YTickMode = 'auto';
    % mark different physical lengths
    if 0
    T1 = axtop.XTick(3)+irf_plot_start_epoch; mean(tint)-diff(tint)*0.2;
    T0 = axtop.XTick(3)+irf_plot_start_epoch; mean(tint);
    T1 = axtop.XTick(1)+irf_plot_start_epoch;
    T_di = di_loc/velocity; t_di = [T1+[0 T_di]' hca.YLim(2)*0.8*[1 1]']; th = text(t_di(2,1)-irf_plot_start_epoch,t_di(2,2)/yratio*0.85,' d_i','color',lengthcolors(1,:),'horizontalalignment','left');
    T_ri = ri_loc/velocity; t_ri = [T1+[0 T_ri]',hca.YLim(2)*0.9*[1 1]']; th = text(t_ri(2,1)-irf_plot_start_epoch,t_ri(2,2)/yratio*0.9,' r_i','color',lengthcolors(2,:),'horizontalalignment','left');
    T_reri = sqrt(ri_loc*re_loc)/velocity; t_reri = [T1+[0 T_reri]',hca.YLim(2)*0.7*[1 1]']; th = text(t_reri(2,1)-irf_plot_start_epoch+0,t_reri(2,2)/yratio*0.75,' (r_er_i)^{1/2}','color',lengthcolors(3,:),'horizontalalignment','left' );
    %T_re = re_loc/velocity; t_re = [T1+[0 T_re]',hca.YLim(2)*0.6*[1 1]']; th = text(t_re(2,1)-irf_plot_start_epoch+0,t_re(2,2)/yratio,'\rho_e','color',lengthcolors(4,:));
    
    hold(hca,'on')
    
    irf_plot(hca,t_di,'linewidth',2,'color',lengthcolors(1,:))
    irf_plot(hca,t_ri,'linewidth',2,'color',lengthcolors(2,:))
    irf_plot(hca,t_reri,'linewidth',2,'color',lengthcolors(3,:))
    %irf_plot(hca,t_re,'linewidth',2)
    hold(hca,'off')  
    end
    grid(hca,'off')    
    hca.XTickLabel = [];
    hca.XLabel.String = '';  
    irf_legend(hca,'d)',[0.02 0.95]);
    velocity;
    
    % Compare Ek and Bz amplitudes
    maxEk = max(dEk(:,2));
    maxEn = max(dEn(:,2));
    maxBz = max(Bz(:,2));
    disp(['flim = ' num2str(flim) '  maxEk = ' num2str(maxEk) ' mv/m,  maxEn = ' num2str(maxEn) ' mv/m,  maxBz = ' num2str(maxBz) ' nT,  maxEk/maxBz = ' num2str(maxEk*1e-3/maxBz/1e-9*1e-3,'%.0f') ' km/s,  maxEk/maxBz/c = ' num2str(maxEk*1e-3/maxBz/1e-9/units.c,'%.3f')])
    
    % high f phi
    flim = 0.5;
    tool.single_event   
    hca = irf_panel('Phi: high f');
    hca.ColorOrder = colors_phi;
    % move panel slightly downwards to make room for km labels on top
    hca.Position = hca.Position + [0 -0.02 0 0];
    irf_plot(hca,{phi_E(:,[1 1+i_v]),phi_B},'comp')
    hca.YLabel.String = '\Phi [V]'; 
    hleg = irf_legend(hca,phi_labels,[0.99 0.95]);
    %hleg = irf_legend(hca,{'\Phi_{\delta E{\perp}}','\Phi_{\delta B{||}}'},[0.99 0.95]);
    %hleg = irf_legend(hca,{'\Phi_{\delta E}','\Phi_{\delta B}'},[0.99 0.95]);
    hleg(1).Color = colors_phi(1,:);
    hleg(2).Color = colors_phi(2,:);
    %irf_legend(hca,{['v_{ph}=' num2str(velocity,'%.0f') ' km/s']},[0.02 0.95])
    irf_legend(hca,{['\omega > ' num2str(flim,'%.2f') ' \Omega_{LH}']},[0.02 0.1])
    irf_zoom(hca,'x',tint)
    irf_zoom(hca,'y',[-15 15]) 
    axtop = add_length_on_top(hca,velocity,0.25);    
    % add B labels on right
    yratio = phi_B(100,2)/Bz(100,2);            
    axtop.YLim = hca.YLim/yratio;
    axtop.YLabel.String = '\delta B_{||} [nT]'; 
    axtop.YTickLabelMode = 'auto';
    axtop.YTickMode = 'auto';
    % mark different physical lengths
    hca = irf_panel('Phi: high f');
    %T1 = axtop.XTick(4)+irf_plot_start_epoch; mean(tint)-diff(tint)*0.2;
    %T0 = axtop.XTick(4)+irf_plot_start_epoch; mean(tint);
    %2
    %T_di = di_loc/velocity; t_di = [T1+[0 T_di]' hca.YLim(2)*0.8*[1 1]']; %th = text(t_di(2,1)-irf_plot_start_epoch,t_di(2,2)/yratio,'d_i')
    %T_re = re_loc/velocity; t_re = [T1+[0 T_re]',hca.YLim(2)*0.9*[1 1]'];
    %T_reri = sqrt(ri_loc*re_loc)/velocity; t_reri = [T1+[0 T_reri]',hca.YLim(2)*0.7*[1 1]'];
    %T_ri = ri_loc/velocity; t_ri = [T1+[0 T_ri]',hca.YLim(2)*0.9*[1 1]'];
    %hold(hca,'on')
    %irf_plot(hca,t_ri,'linewidth',2)
    %irf_plot(hca,t_di,'linewidth',2)
    %irf_plot(hca,t_reri,'linewidth',2)
    %irf_plot(hca,t_re,'linewidth',2,'color',lengthcolors(4,:)); t_re = [T1+[0 T_re]',hca.YLim(2)*0.9*[1 1]']; th = text(t_re(2,1)-irf_plot_start_epoch+0,t_re(2,2)/yratio*0.92,' r_e','color',lengthcolors(4,:));
    %hold(hca,'off')    
    grid(hca,'off')
    irf_legend(hca,'e)',[0.02 0.95]);
    velocity;
    
    irf_plot_axis_align
    irf_plot_zoomin_lines_between_panels(h(3),h(5))
    hca = irf_panel('delete for space'); delete(hca) 
    
    % Compare Ek and Bz amplitudes
     maxEk = max(dEk(:,2));
    maxEn = max(dEn(:,2));
    maxBz = max(Bz(:,2));
    disp(['flim = ' num2str(flim) '  maxEk = ' num2str(maxEk) ' mv/m,  maxEn = ' num2str(maxEn) ' mv/m,  maxBz = ' num2str(maxBz) ' nT,  maxEk/maxBz = ' num2str(maxEk*1e-3/maxBz/1e-9*1e-3,'%.0f') ' km/s,  maxEk/maxBz/c = ' num2str(maxEk*1e-3/maxBz/1e-9/units.c,'%.3f')])
end
%% Spectra figure
if 0 % Spectra figure for paper, old
    %% Some input
    flim = 0.05; % for correlation and velocity
    ffilt = 0.01*flh_loc; % for plotting spectra
    plotFilt = 0;
    doSmooth = 1;
    nantreatment = 'interp';
    v_factor = 1; % multiply v with this to adjust for other time intervals than tint
    fn = 200; % number of f intervals for smoothing
    
    tool.single_event  
    disp(['C = ' num2str(corr_dir(i_dir)) ',  v = ' num2str(velocity,'%.0f') ' x [' num2str(direction) '] km/s'])
    
    filterorder = 3;
    tint_fft = tint;%_ov;
    tint_str= [datestr(irf_time(tint_fft(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint_fft(2),'epoch>datenum'),'hh:mm:ss.fff')];
    
    % Make spectra    
    locE = irf_tlim(gsmE3,tint_fft);
    locEk = locE; locE(:,2) = locE(:,2:4)*direction'; % includes NaNs
    % find nans
    EisNaN = isnan(locE(:,2));
    nanstart = find(EisNaN>0,1,'first');
    nanstop = find(EisNaN>0,1,'last');
    nanall = find(EisNaN>0);  
    
    switch nantreatment
        case 'interp'
        case 'zero'
        case 'cut'
    end
    Ekk = irf_dot(gsmE3(~isnan(gsmE3(:,2)),:),direction); % Take away NaNs so spectra can be made for longer time interval    
    maxEkk = max(Ekk(:,1));
    Bzz = irf_dot(gsmB3,z(i_dir,:));    
    maxBzz = max(Bzz(:,1));
    
    locEk = irf_tlim(Ekk,tint_fft);
    locEkFilt = irf_filt(locEk,ffilt,0,450,filterorder);
    locIntEdt = irf_integrate(locEk);  
    locIntEdtFilt = irf_integrate(locEkFilt);  
    
    locBz = irf_tlim(Bzz,tint_fft);
    locBzFilt = irf_filt(locBz,ffilt,0,450,filterorder);
    locIntBdtFilt = irf_integrate(locBzFilt);  
    locIntBdt = irf_integrate(locBz);  
        
    pfftIntEdt = irf_powerfft(locIntEdt,size(locIntEdt,1),450,0.5); % 
    pfftBz = irf_powerfft(locBz,size(locBz,1),450,0.5);
    
    pfftIntEdtFilt = irf_powerfft(locIntEdtFilt,size(locIntEdtFilt,1),450,0.5); % 
    pfftBzFilt = irf_powerfft(locBzFilt,size(locBzFilt,1),450,0.5);
    pfftIntBdtFilt = irf_powerfft(locIntBdtFilt,size(locBzFilt,1),450,0.5);
    pfftIntBdt = irf_powerfft(locIntBdt,size(locBz,1),450,0.5);
    
    pfftIntBdtFilt = cn_smooth_fft(pfftIntBdtFilt,100);
    pfftIntBdt = cn_smooth_fft(pfftIntBdt,100);
    pfftIntEdtFilt = cn_smooth_fft(pfftIntEdtFilt,100);
    pfftIntEdt = cn_smooth_fft(pfftIntEdt,100);
    if 0
    loglog(pfftIntBdtFilt.f,pfftIntBdtFilt.p{1},pfftIntBdt.f,pfftIntBdt.p{1},...
           pfftIntEdtFilt.f,pfftIntEdtFilt.p{1},pfftIntEdt.f,pfftIntEdt.p{1},...
           pfftEkFilt.f,pfftEkFilt.p{1},pfftIntEdt.f,pfftIntEdt.p{1});
    legend('Bfilt','B','Efilt','E')
    end
    %%
    waveIntEdt = irf_wavelet(locEk,'wavelet_width',5*5.36);
    if 0 % plot wavelet
        %%
        figure(76)
        waveIntEdt = irf_wavelet(locEk,'wavelet_width',5*5.36);
        pcolor(waveIntEdt.t,waveIntEdt.f,log10(waveIntEdt.p{1}'))
        hca=gca;
        shading flat;
        hca.YScale = 'log';
    end        
    
    phiBscale = B0*1e-9/(n_loc*1e6*e*4*pi*1e-7);
    phiEscale = velocity*v_factor;
    
    %lambda = velocity./pfftEk.f; % km    
    %wavenumber = 2*pi./lambda; % 1/km
    %wavenumber = pfftEk.f*2*pi/velocity;

    if 0 % make fits to slopes    
        % dBz 
        bBreak = 150; % kdi=150
        bLow = 10;
        yy = pfftBz.p{1};
        xx = pfftBz.f*2*pi/velocity*di_loc;

        yy1 = yy(xx>bLow)'; yy1 = yy1(xx<bBreak);
        xx1 = xx(xx>bLow); xx1 = xx1(xx<bBreak); 
        P1 = polyfit(log(xx1),log(yy1),1);
        yfit1log = P1(1)*xx1+P1(2);
        yfit1 = exp(P1(1)*log(xx1)+P1(2));
        yfit1fun = @(x) exp(P1(1)*log(x)+P1(2));

        yy2 = yy(xx>bBreak)';
        xx2 = xx(xx>bBreak);        
        P2 = polyfit(log(xx2),log(yy2),1);
        yfit2 = exp(P2(1)*log(xx2)+P2(2));
        yfit2fun = @(x) exp(P2(1)*log(x)+P2(2));

        hold on;
        loglog(xx1,yfit1,'k-',xx2,yfit2,'r-');
    end
    fig=figure(40); fig.NextPlot = 'replace';
    plotE = [locEk(:,1) locEk(:,2)*phiEscale];
    plotB = [locBz(:,1) locBz(:,2)*phiBscale];    
    plotE = [locIntEdtFilt(:,1) locIntEdtFilt(:,2)*phiEscale];
    plotB = [locBzFilt(:,1) locBzFilt(:,2)*1e-9*phiBscale];
    
    hca= subplot(1,1,1);
    irf_plot(hca,{plotE,plotB},'comp');
    legend(hca,'\phi_E','\phi_B')
    
    % Plot figure
    figure(41)
    fontsize = 14;
    klims = [10^-1 7*10^2];    
    nrows = 1;
    
    EFWfilter = 180;
    kfilter = EFWfilter*2*pi/velocity*di_loc;
    lengthcolors;
    
    hca = subplot(nrows,1,1);
    
    %hca.ColorOrder = [colors_phi(1:2,:); lengthcolors];
    if plotFilt
        plotEfft = pfftIntEdtFilt;
        plotBfft = pfftBzFilt;
    else
        plotEfft = pfftIntEdt;
        plotBfft = pfftBz;        
    end
    if doSmooth % Make spectra smooth by averaging        
        plotEfft = cn_smooth_fft(plotEfft,fn);
        plotBfft = cn_smooth_fft(plotBfft,fn);
    end
    
    plotBf = plotBfft.f;
    plotEf = plotEfft.f;
    plotB = plotBfft.p{1}*1e-18*phiBscale^2;
    plotE = plotEfft.p{1}*phiEscale^2;
    
    lines = loglog(hca,plotBfft.f*2*pi/velocity*di_loc,plotB,...
                       plotEfft.f*2*pi/velocity*di_loc,plotE);
    lines(1).Color = colors_phi(2,:);
    lines(2).Color = colors_phi(1,:);
    
    
    hold(hca,'on')
    lines = loglog(di_loc/ri_loc*[1 1],hca.YLim,... 
                   di_loc/sqrt(ri_loc*re_loc)*[1 1],hca.YLim,...
                   di_loc/re_loc*[1 1],hca.YLim);
    lines(1).Color = lengthcolors(2,:);
    lines(2).Color = lengthcolors(3,:);
    lines(3).Color = lengthcolors(4,:);
    %lines = loglog(kfilter*[1 1],hca.YLim,'k--');
    hold(hca,'off')
           
    hca.ColorOrder = lengthcolors;
    %hca.XLim = klims;       
    %h_leg = legend(hca,{'B_{||}\times(B_0/n_ee\mu_0)','\int E_k dt\times v'},'location','northeast','box','off','fontsize',fontsize);
    h_leg = legend(hca,{'\Phi_B','\Phi_E'},'location','northeast','box','off','fontsize',fontsize);
    
    hca.YLabel.String = '[V^2 /Hz]';
    hca.XLabel.String = 'k_{\perp}d_i';
    hca.Title.String = tint_str;
    % Add length scale legends directly to figure
    textposition = 10^-6;hca.YLim(1)+diff(log10(hca.YLim))*0.1; 
    th = text(di_loc/ri_loc,textposition,' k_{\perp}r_i=1','color',lengthcolors(2,:),'fontsize',fontsize);
    th = text(di_loc/sqrt(ri_loc*re_loc),textposition,' k_{\perp}(r_er_i)^{1/2}=1','color',lengthcolors(3,:),'fontsize',fontsize);
    th = text(di_loc/re_loc,textposition,' k_{\perp}r_e=1','color',lengthcolors(4,:),'fontsize',fontsize);
    
    if plotFilt, irf_legend(hca,{['f > ' num2str(ffilt,'%.4f') ' = ' num2str(ffilt/flh_loc,'%.3f') ' f_{LH},  elliptic filter order = ' num2str(filterorder)]},[0.02 0.1]); end
    hca.FontSize = fontsize;
    hca.XLim = [1e-1 1e3];
    hca.YLim = [1e-10 1e8];
end
if 0 % Spectra figure for paper
    %% Some input
    flim = 0.5; % for correlation and velocity
    
    tint_ov = toepoch([2007 09 02 15 46 30; 2007 09 02 15 48 50])';
    tint_ov1 = toepoch([2007 09 02 15 46 30; 2007 09 02 15 47 36])';
    tint_ov2 = toepoch([2007 09 02 15 47 41; 2007 09 02 15 48 50])';
    tint_ov_thinner = toepoch([2007 09 02 15 47 05; 2007 09 02 15 48 04])';
    tint_ov_thinner2 = toepoch([2007 09 02 15 47 05; 2007 09 02 15 48 10])';
    tint = toepoch([2007 09 02 15 47 45.3; 2007 09 02 15 47 46.7])';        
    sc = 3;
    csys = 'gsm';
    
    tool.single_event  
    disp(['C = ' num2str(corr_dir(i_dir)) ',  v = ' num2str(velocity,'%.0f') ' x [' num2str(direction) '] km/s'])
    phiBscale = B0*1e-9/(n_loc*1e6*e*4*pi*1e-7);
    phiEscale = velocity*v_factor;
    
    ffilt = 0.02*flh_loc; % for plotting spectra
    plotFilt = 0;
    doSmooth = 1;
    nantreatment = 'interp';
    v_factor = 1; % multiply v with this to adjust for other time intervals than tint
    fn = 100; % number of f intervals for smoothing    
    filterorder = 3;
    tint_fft = tint_ov;    
    %tint_str= [datestr(irf_time(tint_fft(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint_fft(2),'epoch>datenum'),'hh:mm:ss.fff')];
    tint_str= [datestr(irf_time(tint_fft(1),'epoch>date'),'HH:MM:SS') ' - ' datestr(irf_time(tint_fft(2),'epoch>date'),'HH:MM:SS')];
    
    % Make spectra    
    locE = irf_tlim(gsmE3,tint_fft);
    locEk = locE(:,1:2); locEk(:,2) = locE(:,2:4)*direction'; % includes NaNs
    
    % find nans
    EisNaN = isnan(locE(:,2));
    nanstart = find(EisNaN>0,1,'first');
    nanstop = find(EisNaN>0,1,'last');
    nanall = find(EisNaN>0);  
    
    % replace nans
    switch nantreatment
        case 'interp' % do linear interpolation in nan interval 
            locEkNanInterp = locEk;    
            E1 = locEkNanInterp(nanstart-1,2);
            E2 = locEkNanInterp(nanstop+1,2);    
            Einterp = E1+(nanall-nanall(1))/(nanall(end)-nanall(1))*(E2-E1); % linear interpolation
            locEkNanInterp = locEk; locEkNanInterp(nanall,2) = Einterp;    
            pfftEkNanInterp = irf_powerfft(locEkNanInterp,size(locEkNanInterp,1),450,0.5);
            locIntEkNanInterp = irf_integrate(locEkNanInterp);  
            pfftIntEkNanInterp = irf_powerfft(locIntEkNanInterp,size(locIntEkNanInterp,1),450,0.5);
            locIntEkNanInterpFilt = irf_filt(locIntEkNanInterp,ffilt,0,450,filterorder);
            pfftIntEkNanInterpFilt = irf_powerfft(locIntEkNanInterpFilt,size(locIntEkNanInterpFilt,1),450,0.5);    
            locIntEdt = locIntEkNanInterp;
            locIntEdtFilt = locIntEkNanInterpFilt;
            pfftIntEdt = pfftIntEkNanInterp;
            pfftIntEdtFilt = pfftIntEkNanInterpFilt;
        case 'zero' % set all nans to zero 
            locEkNanZero = locEk(:,[1 2]);
            locEkNanZero(nanall,2) = locEkNanZero(nanall,2)*0;
            pfftEkNanZero = irf_powerfft(locEkNanZero,size(locEkNanZero,1),450,0.5);
            locIntEkNanZero = irf_integrate(locEkNanZero);  
            pfftIntEkNanZero = irf_powerfft(locIntEkNanZero,size(locIntEkNanZero,1),450,0.5);    
            locIntEkNanZeroFilt = irf_filt(locIntEkNanZero,ffilt,0,450,filterorder);
            pfftIntEkNanZeroFilt = irf_powerfft(locIntEkNanZeroFilt,size(locIntEkNanZeroFilt,1),450,0.5);   
            locIntEdt = locIntEkNanZero;
            locIntEdtFilt = locIntEkNanZeroFilt;
            pfftIntEdt = pfftIntEkNanZero;
            pfftIntEdtFilt = pfftIntEkNanZeroFilt;
        case 'cut' % cut off from time series 
            locEkNanCut = locEk; locEkNanCut(nanall,:) = [];
            pfftEkNanCut = irf_powerfft(locEkNanCut,size(locEkNanCut,1),450,0.5);
            locIntEkNanCut = irf_integrate(locEkNanCut);  
            pfftIntEkNanCut = irf_powerfft(locIntEkNanCut,size(locIntEkNanCut,1),450,0.5);
            locIntEkNanCutFilt = irf_filt(locIntEkNanCut,ffilt,0,450,filterorder);
            pfftIntEkNanCutFilt = irf_powerfft(locIntEkNanCutFilt,size(locIntEkNanCutFilt,1),450,0.5);    
            locIntEdt = locIntEkNanCut;
            locIntEdtFilt = locIntEkNanCutFilt;
            pfftIntEdt = pfftIntEkNanCut;
            pfftIntEdtFilt = pfftIntEkNanCutFilt;
    end
    
    Bzz = irf_dot(gsmB3,z(i_dir,:));                   
    locBz = irf_tlim(Bzz,tint_fft);
    locBzFilt = irf_filt(locBz,ffilt,0,450,filterorder);           
    pfftBz = irf_powerfft(locBz,size(locBz,1),450,0.5);
    pfftBzFilt = irf_powerfft(locBzFilt,size(locBz,1),450,0.5);
    
    % Smoothen spectra
    if doSmooth
        clear pffts        
        pffts = who('pfft*');
        for ii = 1:numel(pffts)
            disp(['Smoothing ' pffts{ii}])
            eval([pffts{ii} ' = cn_smooth_fft(' pffts{ii} ',fn);'])
        end
    end
    %pfftIntBdtFilt = cn_smooth_fft(pfftIntBdtFilt,100);
    %pfftIntBdt = cn_smooth_fft(pfftIntBdt,100);
    %pfftIntEdtFilt = cn_smooth_fft(pfftIntEdtFilt,100);
    %pfftIntEdt = cn_smooth_fft(pfftIntEdt,100);  
    
    % Plot wavelet of field, for inspection    
    if 0 % plot wavelet
        %%
        waveE = irf_wavelet(locEk,'wavelet_width',1*5.36);
        waveB = irf_wavelet(locBz,'wavelet_width',1*5.36);
        waveIntEdt = irf_wavelet(locIntEdt,'wavelet_width',1*5.36);
        %%
        figure(76) 
        hca = subplot(3,1,1);
        pcolor(hca,waveE.t,waveE.f,log10(waveE.p{1}'))
        hc = colorbar;
        shading flat;
        hca.YScale = 'log';
        hc.YLabel.String = 'Power [(mV/m)^2/Hz]';
        hca.YLabel.String = 'f [Hz]';
        hca.Title.String = ['Electric field'];
        irf_timeaxis(hca)
        
        hca = subplot(3,1,2);
        pcolor(hca,waveIntEdt.t,waveIntEdt.f,log10(waveIntEdt.p{1}'))
        hc = colorbar;
        shading flat;
        hca.YScale = 'log';
        hc.YLabel.String = 'Power [(s*mV/m)^2/Hz]';
        hca.YLabel.String = 'f [Hz]';
        hca.Title.String = ['Integrated electric field, nantreatment = ' nantreatment];
        irf_timeaxis(hca)
        
        hca = subplot(3,1,3);
        pcolor(hca,waveB.t,waveB.f,log10(waveB.p{1}'))
        hc = colorbar;
        shading flat;
        hca.YScale = 'log';
        hc.YLabel.String = 'Power [(nT)^2/Hz]';
        hca.YLabel.String = 'f [Hz]';
        hca.Title.String = ['Magnetic field'];
        irf_timeaxis(hca)
    end   
         
    % Plot the integrated fields for comparison %%%%%%%%%%%%%%
    fig=figure(42); %fig.NextPlot = 'replace';
    plotIntE = [locIntEdt(:,1) locIntEdt(:,2)*phiEscale];
    plotE = [locEk(:,1) locEk(:,2)];
    plotB = [locBz(:,1) locBz(:,2)*1e-9*phiBscale];    
    plotEFilt = [locIntEdtFilt(:,1) locIntEdtFilt(:,2)*phiEscale];
    plotBFilt = [locBzFilt(:,1) locBzFilt(:,2)*1e-9*phiBscale];
    
    hca = subplot(3,1,1); h(1) = hca;
    irf_plot(hca,{plotIntE,plotB},'comp');    
    legend(hca,'\phi_E','\phi_B')
    irf_pl_mark(hca,[locEk(nanstart,1) locEk(nanstop,1)]);
    irf_pl_mark(hca,tint,'red');
    hca.YLabel.String = '\phi [V]';
    hca.Title.String = ['f_v>' num2str(flim) 'f_{LH},  v = ' num2str(velocity,'%.0f') ' km/s'];
    irf_zoom(hca,'x',tint_fft)
    
    hca= subplot(3,1,2);  h(2) = hca;
    irf_plot(hca,{plotEFilt,plotBFilt},'comp');    
    legend(hca,'\phi_E','\phi_B')
    irf_pl_mark(hca,[locEk(nanstart,1) locEk(nanstop,1)]);
    irf_pl_mark(hca,tint,'red');
    hca.YLabel.String = ['\phi_{f> ' num2str(ffilt) '} [V]']; 
    irf_zoom(hca,'x',tint_fft)    
    
    hca = subplot(3,1,3); h(3) = hca;
    irf_plot(hca,plotE);    
    legend(hca,'E_k')
    irf_pl_mark(hca,[locEk(nanstart,1) locEk(nanstop,1)]);
    irf_pl_mark(hca,tint,'red');
    hca.YLabel.String = 'E [mV/m]';
    hca.Title.String = ['f_v>' num2str(flim) 'f_{LH},  v = ' num2str(velocity,'%.0f') ' km/s'];
    irf_zoom(hca,'x',tint_fft)
    
    linkaxes(h,'x');
        
    % Plot figure with spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fig = figure(43);
    fig.Position = [589 450 538 248];
    fontsize = 14;
    klims = [10^-1 7*10^2];    
    nrows = 1;
    
    EFWfilter = 180;
    kfilter = EFWfilter*2*pi/velocity*di_loc;
    lengthcolors;
    
    hca = subplot(nrows,1,1);
    
    %hca.ColorOrder = [colors_phi(1:2,:); lengthcolors];
    if plotFilt
        plotEfft = pfftIntEdtFilt;
        plotBfft = pfftBzFilt;
    else
        plotEfft = pfftIntEdt;
        plotBfft = pfftBz;        
    end
    %if doSmooth % Make spectra smooth by averaging        
    %    plotEfft = cn_smooth_fft(plotEfft,fn);
    %    plotBfft = cn_smooth_fft(plotBfft,fn);
    %end
    
    plotBf = plotBfft.f;
    plotEf = plotEfft.f;
    plotB = plotBfft.p{1}*1e-18*phiBscale^2;
    plotE = plotEfft.p{1}*phiEscale^2;
    idxNan = isnan(plotE); 
    plotE(idxNan) = [];
    plotB(idxNan) = [];
    plotEf(idxNan) = [];
    plotBf(idxNan) = [];
    
    lines = loglog(hca,plotBf*2*pi/velocity*di_loc,plotB,...
                       plotEf*2*pi/velocity*di_loc,plotE);
    lines(1).Color = colors_phi(2,:);
    lines(2).Color = colors_phi(1,:);
    
    hca.YLabel.String = 'Power [V^2 /Hz]';
    hca.XLabel.String = 'k_{\perp}d_i';
    hca.Title.String = tint_str;
    hca.XLim = [2e-2 1e3];
    hca.YLim = [1e-5 1e9];
    tickstep = 1;
    %hca.YTick = 10.^(log10(hca.YLim(1)):tickstep:log10(hca.YLim(2)));
    hca.YTick = 10.^(-16:2:16);
    
    hold(hca,'on')
    lines = loglog(di_loc/ri_loc*[1 1],hca.YLim,... 
                   di_loc/sqrt(ri_loc*re_loc)*[1 1],hca.YLim,...
                   di_loc/re_loc*[1 1],hca.YLim);
    lines(1).Color = lengthcolors(2,:);
    lines(2).Color = lengthcolors(3,:);
    lines(3).Color = lengthcolors(4,:);
    %lines = loglog(kfilter*[1 1],hca.YLim,'k--');
    hold(hca,'off')
           
    hca.ColorOrder = lengthcolors;
    %hca.XLim = klims;       
    %h_leg = legend(hca,{'B_{||}\times(B_0/n_ee\mu_0)','\int E_k dt\times v'},'location','northeast','box','off','fontsize',fontsize);
    h_leg = legend(hca,{'\Phi_B','\Phi_E'},'location','northeast','box','off','fontsize',fontsize);
    
   
    
    % Add length scale legends directly to figure
    halign = 'left';
    textposition = 3*10^-4;hca.YLim(1)+diff(log10(hca.YLim))*0.1; 
    th = text(di_loc/ri_loc,textposition,' k_{\perp}r_i=1','color',lengthcolors(2,:),'fontsize',fontsize,'horizontalalignment',halign);
    th = text(di_loc/sqrt(ri_loc*re_loc),textposition,' k_{\perp}(r_er_i)^{1/2}=1','color',lengthcolors(3,:),'fontsize',fontsize,'horizontalalignment',halign);
    th = text(di_loc/re_loc,textposition,' k_{\perp}r_e=1','color',lengthcolors(4,:),'fontsize',fontsize,'horizontalalignment',halign);
    
    if plotFilt, irf_legend(hca,{['f > ' num2str(ffilt,'%.4f') ' = ' num2str(ffilt/flh_loc,'%.3f') ' f_{LH},  elliptic filter order = ' num2str(filterorder)]},[0.02 0.1]); end
    hca.FontSize = fontsize;
    
end
if 1 % Spectra figure for paper with integration in Fourier space
    %% Some input
    flim = 0.1; % for correlation and velocity
    
    tint_ov = toepoch([2007 09 02 15 46 30; 2007 09 02 15 48 50])';
    tint_ov1 = toepoch([2007 09 02 15 46 30; 2007 09 02 15 47 36])';
    tint_ov2 = toepoch([2007 09 02 15 47 41; 2007 09 02 15 48 50])';
    tint_ov_thinner = toepoch([2007 09 02 15 47 05; 2007 09 02 15 48 04])';
    tint_ov_thinner2 = toepoch([2007 09 02 15 47 05; 2007 09 02 15 48 10])';
    tint = toepoch([2007 09 02 15 47 45.3; 2007 09 02 15 47 46.7])';  
    tint_noactivity = toepoch([2007 09 02 15 46 30; 2007 09 02 15 47 00])';            
    sc = 3;
    csys = 'gsm';
    
    tool.single_event  
    disp(['C = ' num2str(corr_dir(i_dir)) ',  v = ' num2str(velocity,'%.0f') ' x [' num2str(direction) '] km/s'])
    
    
    doEn = 0;
    doE45 = 0;
    ffilt = 0.02*flh_loc; % for plotting spectra
    lowpassB = 0.1;
    plotFilt = 0;
    doSmooth = 1;
    nantreatment = 'interp';
    v_factor = 1;320/430; % multiply v with this to adjust for other time intervals than tint
    fn = 150; % number of f intervals for smoothing    
    filterorder = 5;
    tint_fft = tint_ov_thinner2;    
    %tint_str= [datestr(irf_time(tint_fft(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint_fft(2),'epoch>datenum'),'hh:mm:ss.fff')];
    tint_str= [datestr(irf_time(tint_fft(1),'epoch>date'),'HH:MM:SS') ' - ' datestr(irf_time(tint_fft(2),'epoch>date'),'HH:MM:SS')];
    
    phiBscale = B0*1e-9/(n_loc*1e6*e*4*pi*1e-7);
    phiEscale = velocity*v_factor;
    phiEscale = 1; % omega/v = 2*pi*f/v
    
    phi_labels_all = {{'\Phi_{\delta B_{||}}','\Phi_{\delta E_{k}}'},...
                      {'\Phi_{\delta B{||}}','\Phi_{\delta E{k}}'},...
                      {'\Phi_{\delta B_{||}}','\Phi_{\delta E_{\perp}}'},...
                      {'\Phi_{\delta B{||}}','\Phi_{\delta E{\perp}}'},...
                      {'\Phi_{\delta B}','\Phi_{\delta E}'},...
                      {'\Phi_{B}','\Phi_{E}'},...
                      {'\Phi_{B_{||}}','\Phi_{E_{k}}'},...
                      {'\Phi_{B{||}}','\Phi_{E{k}}'},...
                      {'\Phi_{B_{||}}','\Phi_{E_{\perp}}'},...
                      {'\Phi_{B{||}}','\Phi_{E{\perp}}'},...
                      {'\Phi_{B_{||}}','\Phi_{E_{n}}'},...
                      {'\Phi_{B_{||}}','\Phi_{E_{45}}'}};
    phi_labels = phi_labels_all{9};
    % Make spectra    
    locE = irf_tlim(gsmE3,tint_fft);    
    if doEn % compare difference if instead integrating En, includes NaNs        
        locEk = locE(:,1:2); locEk(:,2) = locE(:,2:4)*y(i_dir,:)';direction';
        phi_labels = phi_labels_all{11};
    elseif doE45 
        locEk = locE(:,1:2); locEk(:,2) = locE(:,2:4)*direction';
        locEn = locE(:,1:2); locEn(:,2) = locE(:,2:4)*y(i_dir,:)';
        
        %locE45 = [locE(:,1) (locEk(:,2)+locEn(:,2))./sqrt(locEk(:,2).^2+locEn(:,2).^2)];
        dir45 = (direction+y(i_dir,:))/sqrt(2);       
        locE45 = locE(:,1:2); locE45(:,2) = locE(:,2:4)*dir45';
        locEk = locE45;       
        phi_labels = {'\Phi_{B_{||}}','\Phi_{E_{45}}'};
    else        
        locEk = locE(:,1:2); locEk(:,2) = locE(:,2:4)*direction'; % includes NaNs
    end
    locEkFilt = irf_filt(locEk,ffilt,0,450,filterorder); % includes NaNs    
    % find nans
    EisNaN = isnan(locE(:,2));
    nanstart = find(EisNaN>0,1,'first');
    nanstop = find(EisNaN>0,1,'last');
    nanall = find(EisNaN>0); 
    
    
    
    % replace nans, but not nned to integrate!
    if ~isempty(nanall)        
        switch nantreatment
            case 'interp' % do linear interpolation in nan interval 
                locEkNanInterp = locEk;    
                E1 = locEkNanInterp(nanstart-1,2);
                E2 = locEkNanInterp(nanstop+1,2);    
                Einterp = E1+(nanall-nanall(1))/(nanall(end)-nanall(1))*(E2-E1); % linear interpolation
                locEkNanInterp = locEk; locEkNanInterp(nanall,2) = Einterp;    
                pfftEkNanInterp = irf_powerfft(locEkNanInterp,size(locEkNanInterp,1),450,0.5);           
                pfftEk = pfftEkNanInterp;            
            case 'zero' % set all nans to zero 
                locEkNanZero = locEk(:,[1 2]);
                locEkNanZero(nanall,2) = locEkNanZero(nanall,2)*0;
                pfftEkNanZero = irf_powerfft(locEkNanZero,size(locEkNanZero,1),450,0.5);            
                pfftEk = pfftEkNanZero;
            case 'cut' % cut off from time series 
                locEkNanCut = locEk; locEkNanCut(nanall,:) = [];
                pfftEkNanCut = irf_powerfft(locEkNanCut,size(locEkNanCut,1),450,0.5); 
                pfftEk = pfftEkNanCut;
        end
    else
        pfftEk = irf_powerfft(locEk,size(locEk,1),450,0.5); 
    end
    
    pfftPhiE = pfftEk;
    pfftPhiE.p{1} = pfftEk.p{1}.*(velocity*v_factor/2/pi./pfftEk.f').^2;
    
    %Bzz = irf_dot(gsmB3,z(i_dir,:));                   
    Bzz = irf_abs(gsmB3);
    locBz = irf_tlim(Bzz(:,[1 5]),tint_fft);
    BDC = irf_filt(locBz,0,lowpassB,450,3);
    locPhiScale = irf_multiply(1e-9/(1e6*e*4*pi*1e-7),BDC,1,peaNe3hf,-1);
    locPhiB = irf_multiply(1e-9,locBz,1,locPhiScale,1); % V
    locBzFilt = irf_filt(locBz,ffilt,0,450,filterorder);           
    pfftBz = irf_powerfft(locBz,size(locBz,1),450,0.5);
    pfftPhiBbetter = irf_powerfft(locPhiB,size(locPhiB,1),450,0.5);
    
    figure(44);irf_plot({locBz,BDC})
    
    %pfftBzFilt = irf_powerfft(locBzFilt,size(locBz,1),450,0.5);
    %pfftEk = pfftEkTs; pfftEk.p{1} = irf.nanmean(pfftEkTs.p{1},1); pfftEk.t = irf.nanmean(pfftEkTs.t,1);
    %nfft = 2048;
    
    
    
    
    if 0 % use local values for phiB
        pfftPhiB = pfftBz; 
        phiBscale = B0*1e-9/(n_loc*1e6*e*4*pi*1e-7);
        pfftPhiB.p{1} = pfftPhiB.p{1}*1e-18*phiBscale^2;
    elseif 0% use a single set of values for B0 and n
        pfftPhiB = pfftBz; 
        phiBscale = irf_multiply(1e-9/(1e6*e*4*pi*1e-7),Bzz(:,[1 5]),1,peaNe3hf,-1);
        pfftPhiB.p{1} = pfftPhiB.p{1}*1e-18*phiBscale^2;
        %%
        if 0 % plot relative difference
            irf_plot({phiBscale,[phiBscale(:,1) phiBscale(:,2)/(B0*1e-9/(n_loc*1e6*e*4*pi*1e-7))-1]})
            irf_pl_mark(gca,tint_fft)
        end
    else
        pfftPhiB = pfftPhiBbetter; 
    end
    
    % Smoothen spectra
    if doSmooth
        clear pffts        
        pffts = who('pfft*');
        for ii = 1:numel(pffts)
            disp(['Smoothing ' pffts{ii}])
            eval([pffts{ii} ' = cn_smooth_fft(' pffts{ii} ',fn);'])
        end
    end    
    
    % Plot wavelet of field, for inspection    
    if 0 % plot wavelet
        %%
        waveE = irf_wavelet(locEk,'wavelet_width',1*5.36);
        waveB = irf_wavelet(locBz,'wavelet_width',1*5.36);
        waveIntEdt = irf_wavelet(locIntEdt,'wavelet_width',1*5.36);
        %%
        figure(76) 
        hca = subplot(3,1,1);
        pcolor(hca,waveE.t,waveE.f,log10(waveE.p{1}'))
        hc = colorbar;
        shading flat;
        hca.YScale = 'log';
        hc.YLabel.String = 'Power [(mV/m)^2/Hz]';
        hca.YLabel.String = 'f [Hz]';
        hca.Title.String = ['Electric field'];
        irf_timeaxis(hca)
        
        hca = subplot(3,1,2);
        pcolor(hca,waveIntEdt.t,waveIntEdt.f,log10(waveIntEdt.p{1}'))
        hc = colorbar;
        shading flat;
        hca.YScale = 'log';
        hc.YLabel.String = 'Power [(s*mV/m)^2/Hz]';
        hca.YLabel.String = 'f [Hz]';
        hca.Title.String = ['Integrated electric field, nantreatment = ' nantreatment];
        irf_timeaxis(hca)
        
        hca = subplot(3,1,3);
        pcolor(hca,waveB.t,waveB.f,log10(waveB.p{1}'))
        hc = colorbar;
        shading flat;
        hca.YScale = 'log';
        hc.YLabel.String = 'Power [(nT)^2/Hz]';
        hca.YLabel.String = 'f [Hz]';
        hca.Title.String = ['Magnetic field'];
        irf_timeaxis(hca)
    end                
        
    % Plot figure with spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fig = figure(43);
    fig.Position = [589 450 538 248];
    fontsize = 14;
    klims = [10^-1 7*10^2];    
    nrows = 1;
    
    EFWfilter = 180;
    kfilter = EFWfilter*2*pi/velocity*di_loc;
    lengthcolors;
        
    
    %hca.ColorOrder = [colors_phi(1:2,:); lengthcolors];
    if plotFilt
        plotEfft = pfftEkFilt;
        plotBfft = pfftBzFilt;
    else
        plotEfft = pfftEk;
        plotBfft = pfftBz;        
    end
    %if doSmooth % Make spectra smooth by averaging        
    %    plotEfft = cn_smooth_fft(plotEfft,fn);
    %    plotBfft = cn_smooth_fft(plotBfft,fn);
    %end
   
    plotPhiBf = pfftPhiB.f;
    plotPhiEf = pfftPhiE.f;
    plotPhiBp = pfftPhiB.p{1};
    plotPhiEp = pfftPhiE.p{1};
    % remove values when frequency bin has no input (= NaN)
    % E
    idxNan = isnan(plotPhiEp); 
    plotPhiEf(idxNan) = [];
    plotPhiEp(idxNan) = [];
    idxBNan = isnan(plotPhiBp);
    % B
    plotPhiBf(idxBNan) = [];
    plotPhiBp(idxBNan) = [];    
    % remove values above 180 Hz filter
    % E
    idxE180 = find(plotPhiEf>180);
    plotPhiEf(idxE180) = [];
    plotPhiEp(idxE180) = [];
    % B
    idxB180 = find(plotPhiBf>180); 
    plotPhiBf(idxB180) = [];   
    plotPhiBp(idxB180) = [];
    % remove frequencies which is above the magnetic field noise floor, ~80
    % Hz, see figure that can be plottet a few lines below
    idxBnoise = find(plotPhiBf>80);
    plotPhiBp(idxBnoise) = [];
    plotPhiBf(idxBnoise) = [];
    
    % remove frequencies which are below 0.03 Hz 
    lowestf = 0.03;
    idxBlow = find(plotPhiBf<lowestf);
    plotPhiBp(idxBlow) = [];
    plotPhiBf(idxBlow) = [];
    idxElow = find(plotPhiEf<lowestf);
    plotPhiEp(idxElow) = [];
    plotPhiEf(idxElow) = [];
    
    % plot figure
    hca = subplot(1,1,1);
    lines = loglog(hca,plotPhiBf*2*pi/velocity*di_loc,plotPhiBp,'-',...
                       plotPhiEf*2*pi/velocity*di_loc,plotPhiEp,'-');
    lines(1).Color = colors_phi(2,:);
    lines(2).Color = colors_phi(1,:);
    
    hca.YLabel.String = 'Power [V^2 /Hz]';
   %hca.YLabel.String = 'Power Spectral Density [V^2 /Hz]';
    hca.XLabel.String = 'k_{\perp}d_i';
    hca.Title.String = tint_str;
    hca.XLim = [3e-2 1e3];
    hca.YLim = [1e-6 1e10];
    tickstep = 1;
    %hca.YTick = 10.^(log10(hca.YLim(1)):tickstep:log10(hca.YLim(2)));
    hca.YTick = 10.^(-16:2:16);
    
    hold(hca,'on')
    lines = loglog(di_loc/ri_loc*[1 1],hca.YLim,... 
                   di_loc/sqrt(ri_loc*re_loc)*[1 1],hca.YLim,...
                   di_loc/re_loc*[1 1],hca.YLim);
    lines(1).Color = lengthcolors(2,:);
    lines(2).Color = lengthcolors(3,:);
    lines(3).Color = lengthcolors(4,:);
    %lines = loglog(kfilter*[1 1],hca.YLim,'k--');
    hold(hca,'off')
           
    hca.ColorOrder = lengthcolors;
    %hca.XLim = klims;       
    %h_leg = legend(hca,{'B_{||}\times(B_0/n_ee\mu_0)','\int E_k dt\times v'},'location','northeast','box','off','fontsize',fontsize);
    h_leg = legend(hca,phi_labels,'location','northeast','box','off','fontsize',fontsize);
    
   
    
    % Add length scale legends directly to figure
    halign = 'left';
    textposition = 3*10^-4;hca.YLim(1)+diff(log10(hca.YLim))*0.1; 
    th = text(di_loc/ri_loc,textposition,' k_{\perp}r_i=1','color',lengthcolors(2,:),'fontsize',fontsize,'horizontalalignment',halign);
    th = text(di_loc/sqrt(ri_loc*re_loc),textposition,' k_{\perp}(r_er_i)^{1/2}=1','color',lengthcolors(3,:),'fontsize',fontsize,'horizontalalignment',halign);
    th = text(di_loc/re_loc,textposition,' k_{\perp}r_e=1','color',lengthcolors(4,:),'fontsize',fontsize,'horizontalalignment',halign);
    
    if plotFilt, irf_legend(hca,{['f > ' num2str(ffilt,'%.4f') ' = ' num2str(ffilt/flh_loc,'%.3f') ' f_{LH},  elliptic filter order = ' num2str(filterorder)]},[0.02 0.1]); end
    hca.FontSize = fontsize;
    
    if 0 % check to see if B hits noise floor
        %%
        % Need to compare spectra outside when there is 'nothing' happening
        % and where we have waves
        
        tint_noactivity = toepoch([2007 09 02 15 46 30; 2007 09 02 15 47 00])';            
        tint_ov_thinner2 = toepoch([2007 09 02 15 47 05; 2007 09 02 15 48 10])';
        tint = toepoch([2007 09 02 15 47 45.3; 2007 09 02 15 47 46.7])';        
        
        Bzz = irf_dot(gsmB3,z(i_dir,:));                   
        locBz = irf_tlim(Bzz,tint_ov_thinner);
        pfftBz = irf_powerfft(locBz,size(locBz,1),450,0.5);
        locBzNoact = irf_tlim(Bzz,tint_noactivity);
        pfftBzNoact = irf_powerfft(locBzNoact,size(locBzNoact,1),450,0.5);
        
        locE = irf_tlim(gsmE3,tint_noactivity);
        locEk = locE(:,1:2); locEk(:,2) = locE(:,2:4)*direction'; % includes NaNs
    
        locEkNoact = irf_tlim(locEk,tint_noactivity);
        pfftEkNoact = irf_powerfft(locEkNoact,size(locEkNoact,1),450,0.5);
        pfftEkBNoact = pfftEkNoact; pfftEkBNoact.p{1}=pfftEkNoact.p{1}.*(velocity/2/pi./pfftEkNoact.f').^2/(1e-18*phiBscale^2);
        
        pfftBzNoactSmooth = cn_smooth_fft(pfftBzNoact,fn);
        pfftBzSmooth = cn_smooth_fft(pfftBz,fn);        
        pfftEkBNoactSmooth = cn_smooth_fft(pfftEkBNoact,fn);        
        
        figure(54)
        hca=subplot(3,1,1); h(1) = hca;
        irf_plot(hca,{locBzNoact,locBz},'comp')
        hca.YLabel.String = 'B [nT]';
        hca=subplot(3,1,2); h(2) = hca;
        loglog(hca,pfftBzNoact.f,pfftBzNoact.p{1},pfftBz.f,pfftBz.p{1})
        hca.YLabel.String = 'Log Power spectral Density [nT^2/Hz]';
        hca.XLabel.String = 'Frequency [Hz]';
        
        hca=subplot(3,1,3); h(3) = hca;
        loglog(hca,pfftBzNoactSmooth.f,pfftBzNoactSmooth.p{1},'.',pfftBzSmooth.f,pfftBzSmooth.p{1},'.')%,...
                   %pfftEkBNoactSmooth.f,pfftEkBNoactSmooth.p{1},'.')
        hca.YLabel.String = 'Log Power spectral Density [nT^2/Hz]';
        hca.XLabel.String = 'Frequency [Hz]';
        hca.Title.String = 'Averaged over frequencies';
        
        linkaxes(h,'off')
        
        
    end
end
if 1 % Spectra figure for paper with integration in Fourier space and taking average in time instead of in frequencies
    %% Some input
    flim = 0.5; % for correlation and velocity
    
    tint_ov = toepoch([2007 09 02 15 46 30; 2007 09 02 15 48 50])';
    tint_ov1 = toepoch([2007 09 02 15 46 30; 2007 09 02 15 47 36])';
    tint_ov2 = toepoch([2007 09 02 15 47 41; 2007 09 02 15 48 50])';
    tint_ov_thinner = toepoch([2007 09 02 15 47 05; 2007 09 02 15 48 04])';
    tint_ov_thinner2 = toepoch([2007 09 02 15 47 05; 2007 09 02 15 48 10])';
    tint = toepoch([2007 09 02 15 47 45.3; 2007 09 02 15 47 46.7])';  
    tint_noactivity = toepoch([2007 09 02 15 46 30; 2007 09 02 15 47 00])';            
    sc = 3;
    csys = 'gsm';
    
    tool.single_event  
    disp(['C = ' num2str(corr_dir(i_dir)) ',  v = ' num2str(velocity,'%.0f') ' x [' num2str(direction) '] km/s'])
    phiBscale = B0*1e-9/(n_loc*1e6*e*4*pi*1e-7);
    phiEscale = velocity*v_factor;
    phiEscale = 1; % omega/v = 2*pi*f/v
    
    ffilt = 0.02*flh_loc; % for plotting spectra
    plotFilt = 0;
    doSmooth = 1;
    waveletwidth = 1;
    tint_fft = tint_ov_thinner;    
    tint_str= [datestr(irf_time(tint_fft(1),'epoch>date'),'HH:MM:SS') ' - ' datestr(irf_time(tint_fft(2),'epoch>date'),'HH:MM:SS')];
     
    % Make spectra    
    locE = irf_tlim(gsmE3,tint_fft);
    locEk = locE(:,1:2); locEk(:,2) = locE(:,2:4)*direction'; % includes NaNs
    locEkFilt = irf_filt(locEk,ffilt,0,450,filterorder); % includes NaNs
    
    nfft = 2048;
    pfftEkTs = irf_powerfft(locEk,nfft,450,0); 
    waveEkTs = irf_wavelet(locEk,'wavelet_width',waveletwidth*5.36,'f',[0.1 300],'nf',100);
    
    
    % find nans
    EisNaN = isnan(locE(:,2));
    nanstart = find(EisNaN>0,1,'first');
    nanstop = find(EisNaN>0,1,'last');
    nanall = find(EisNaN>0); 

    % do linear interpolation in nan interval 
    locEkNanInterp = locEk;    
    E1 = locEkNanInterp(nanstart-1,2);
    E2 = locEkNanInterp(nanstop+1,2);    
    Einterp = E1+(nanall-nanall(1))/(nanall(end)-nanall(1))*(E2-E1); % linear interpolation
    locEkNanInterp = locEk; locEkNanInterp(nanall,2) = Einterp;    
    pfftEkNanInterp = irf_powerfft(locEkNanInterp,size(locEkNanInterp,1),450,0.5);           
    pfftEkOnesingle = pfftEkNanInterp;  
        
    % make average
    waveEk = waveEkTs; waveEk.p{1} = irf.nanmean(waveEkTs.p{1},1); waveEk.t = irf.nanmean(waveEkTs.t,1);
    pfftEk = pfftEkTs; pfftEk.p{1} = irf.nanmean(pfftEkTs.p{1},1); pfftEk.t = irf.nanmean(pfftEkTs.t,1);
    %
    if 0 % plot pfft and wavelet 
        %% plot
        figure(99)
        h = irf_plot(3);
        
        isub = 1;
        hca=h(isub); isub=isub+1;
        irf_plot(hca,locEk)
        
        hca=h(isub); isub=isub+1;
        irf_spectrogram(hca,pfftEkTs,'log')
        hca.YScale = 'log';
        shading flat;
        
        hca=h(isub); isub=isub+1;
        irf_spectrogram(hca,waveEkTs)
        hca.YScale = 'log';
        shading flat;
        
        irf_zoom(h,'x',tint_fft)
    end
    %
    pfftPhiE = pfftEk; 
    pfftPhiE = pfftEkOnesingle;
    pfftPhiE.p{1} = pfftPhiE.p{1}.*(velocity/2/pi./pfftPhiE.f').^2;
    
    wavePhiE = waveEk;
    wavePhiE.p{1} = wavePhiE.p{1}.*(velocity/2/pi./wavePhiE.f').^2;
    
    
    Bzz = irf_dot(gsmB3,z(i_dir,:));                   
    locBz = irf_tlim(Bzz,tint_fft);
    locBzFilt = irf_filt(locBz,ffilt,0,450,filterorder);           
    pfftBz = irf_powerfft(locBz,size(locBz,1),450,0.5);
    pfftBzFilt = irf_powerfft(locBzFilt,size(locBz,1),450,0.5);
    
    pfftPhiB = pfftBz; pfftPhiB.p{1} = pfftPhiB.p{1}*1e-18*phiBscale^2;
    pfftPhiB
    % Smoothen spectra
    if doSmooth
        clear pffts        
        pffts = who('pfft*');
        for ii = 1:numel(pffts)
            disp(['Smoothing ' pffts{ii}])
            eval([pffts{ii} ' = cn_smooth_fft(' pffts{ii} ',fn);'])
        end
    end          
        
    % Plot figure with spectra %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    fig = figure(43);
    fig.Position = [589 450 538 248];
    fontsize = 14;
    klims = [10^-1 7*10^2];    
    nrows = 1;
    
    EFWfilter = 180;
    kfilter = EFWfilter*2*pi/velocity*di_loc;
    lengthcolors;
        
    
    %hca.ColorOrder = [colors_phi(1:2,:); lengthcolors];
    if plotFilt
        plotEfft = pfftEkFilt;
        plotBfft = pfftBzFilt;
    else
        plotEfft = pfftEk;
        plotBfft = pfftBz;        
    end
    %if doSmooth % Make spectra smooth by averaging        
    %    plotEfft = cn_smooth_fft(plotEfft,fn);
    %    plotBfft = cn_smooth_fft(plotBfft,fn);
    %end
   
    plotPhiBf = pfftPhiB.f;
    plotPhiEf = pfftPhiE.f;
    plotPhiEfWave = wavePhiE.f;
    plotPhiBp = pfftPhiB.p{1};
    plotPhiEp = pfftPhiE.p{1};
    plotPhiEpWave = wavePhiE.p{1};
    idxNan = isnan(plotPhiEp); 
    idxBNan = isnan(plotPhiBp); 
    plotPhiEf(idxNan) = [];
    plotPhiBf(idxBNan) = [];
    plotPhiEp(idxNan) = [];
    plotPhiBp(idxBNan) = [];
    plotPhiEfWave(idxNan) = [];
    plotPhiEpWave(idxNan) = [];
    idxWaveE180 = find(plotPhiEfWave>180);
    idxE180 = find(plotPhiEf>180);
    idxB180 = find(plotPhiBf>180);
    plotPhiEf(idxE180) = [];
    plotPhiBf(idxB180) = [];
    plotPhiEp(idxE180) = [];
    plotPhiBp(idxB180) = [];
    plotPhiEf(idxWaveE180) = [];
    plotPhiEp(idxWaveE180) = [];
    %
    hca = subplot(1,1,1);
    lines = loglog(hca,plotPhiBf*2*pi/velocity*di_loc,plotPhiBp,'-',...
                       plotPhiEf*2*pi/velocity*di_loc,plotPhiEp,'-',...
                       plotPhiEfWave*2*pi/velocity*di_loc,plotPhiEpWave,'-');
    lines(1).Color = colors_phi(2,:);
    lines(2).Color = colors_phi(1,:);
    lines(3).Color = colors_phi(4,:);
    
    hca.YLabel.String = 'Power [V^2 /Hz]';
    hca.XLabel.String = 'k_{\perp}d_i';
    hca.Title.String = tint_str;
    hca.XLim = [1e-1 1e3];
    hca.YLim = [1e-6 1e11];
    tickstep = 1;
    %hca.YTick = 10.^(log10(hca.YLim(1)):tickstep:log10(hca.YLim(2)));
    hca.YTick = 10.^(-16:2:16);
    
    hold(hca,'on')
    lines = loglog(di_loc/ri_loc*[1 1],hca.YLim,... 
                   di_loc/sqrt(ri_loc*re_loc)*[1 1],hca.YLim,...
                   di_loc/re_loc*[1 1],hca.YLim);
    lines(1).Color = lengthcolors(2,:);
    lines(2).Color = lengthcolors(3,:);
    lines(3).Color = lengthcolors(4,:);
    %lines = loglog(kfilter*[1 1],hca.YLim,'k--');
    hold(hca,'off')
           
    hca.ColorOrder = lengthcolors;
    %hca.XLim = klims;       
    %h_leg = legend(hca,{'B_{||}\times(B_0/n_ee\mu_0)','\int E_k dt\times v'},'location','northeast','box','off','fontsize',fontsize);
    h_leg = legend(hca,{'\Phi_B',['\Phi_{E,nfft=' num2str(nfft) '}'],'\Phi_{E,wavelet}'},'location','northeast','box','off','fontsize',fontsize);
    
   
    
    % Add length scale legends directly to figure
    halign = 'left';
    textposition = 3*10^-4;hca.YLim(1)+diff(log10(hca.YLim))*0.1; 
    th = text(di_loc/ri_loc,textposition,' k_{\perp}r_i=1','color',lengthcolors(2,:),'fontsize',fontsize,'horizontalalignment',halign);
    th = text(di_loc/sqrt(ri_loc*re_loc),textposition,' k_{\perp}(r_er_i)^{1/2}=1','color',lengthcolors(3,:),'fontsize',fontsize,'horizontalalignment',halign);
    th = text(di_loc/re_loc,textposition,' k_{\perp}r_e=1','color',lengthcolors(4,:),'fontsize',fontsize,'horizontalalignment',halign);
    
    if plotFilt, irf_legend(hca,{['f > ' num2str(ffilt,'%.4f') ' = ' num2str(ffilt/flh_loc,'%.3f') ' f_{LH},  elliptic filter order = ' num2str(filterorder)]},[0.02 0.1]); end
    hca.FontSize = fontsize;            
end
if 0 % Spectra example figure
    %% Do data
    flim = 0.01;
    tool.single_event          
    Ekk = irf_dot(gsmE3(~isnan(gsmE3(:,2)),:),direction);
    %Ekk = irf_filt(Ekk,0.5,0,450,5);
    Bzz = irf_dot(gsmB3,z(i_dir,:));
    
    tint_fft = tint+2*[-1 10];
    tint_fft = tint_ov;
    locEk = irf_tlim(Ekk,tint_fft);
    locPhik = irf_integrate(locEk); locPhik(:,2) = locPhik(:,2)*velocity;
    locBz = irf_tlim(Bzz,tint_fft);
    
    pfftEk = irf_powerfft(locEk,size(locEk,1),450,0.5);    
    pfftPhik = irf_powerfft(locPhik,size(locPhik,1),450,0.5);
    pfftBz = irf_powerfft(locBz,size(locBz,1),450,0.5);
    EBscale = B0*1e-9/(n_loc*1e6*e*4*pi*1e-7);
    
    pfftPhiE = irf_powerfft(phi_E(:,[1 1+i_v]),size(phi_E(:,[1 1+i_v]),1),450,0.5);
    pfftIntEdt = irf_powerfft(intEdt(:,[1 1+i_dir]),size(intEdt(:,[1 1+i_dir]),1),450,0.5);
    pfftPhiB = irf_powerfft(phi_B,size(phi_B,1),450,0.5);
    
    
    
    lambda = velocity./pfftEk.f; % km    
    wavenumber = 2*pi./lambda; % 1/km
    wavenumber = pfftEk.f*2*pi/velocity;

    %% make fits to slopes
    % dBz 
    bBreak = 150; % kdi=150
    bLow = 10;
    yy = pfftBz.p{1};
    xx = pfftBz.f*2*pi/velocity*di_loc;
    
    yy1 = yy(xx>bLow)'; yy1 = yy1(xx<bBreak);
    xx1 = xx(xx>bLow); xx1 = xx1(xx<bBreak); 
    P1 = polyfit(log(xx1),log(yy1),1);
    yfit1log = P1(1)*xx1+P1(2);
    yfit1 = exp(P1(1)*log(xx1)+P1(2));
    yfit1fun = @(x) exp(P1(1)*log(x)+P1(2));
    
    yy2 = yy(xx>bBreak)';
    xx2 = xx(xx>bBreak);        
    P2 = polyfit(log(xx2),log(yy2),1);
    yfit2 = exp(P2(1)*log(xx2)+P2(2));
    yfit2fun = @(x) exp(P2(1)*log(x)+P2(2));
    
    hold on;
    loglog(xx1,yfit1,'k-',xx2,yfit2,'r-');
    
    %% Plot
    klims = [10^-1 5*10^2];
    %klims = [80*10^-1 5*10^2];
    nrows = 4;
    hca = subplot(nrows,1,1);
    loglog(hca,pfftEk.f/flh_loc,pfftEk.p{1},...
               pfftBz.f/flh_loc,pfftBz.p{1},...
               pfftEk.f/flh_loc,pfftEk.p{1}./(pfftEk.f'.^2))        
    hca.XLim = [10^-2 2*10^1];
    legend(hca,{'E_k','B_{||}','E_k/f'},'location','best')
    hca.YLabel.String = 'Power';
    hca.XLabel.String = 'f/f_LH';
    
    %hca = subplot(nrows,1,2);
    %semilogy(hca,pfftPhik.f/flh_loc,pfftPhik.p{1},...
    %             pfftBz.f/flh_loc,pfftBz.p{1}*1e-18*EBscale^2,...
    %             pfftIntEdt.f/flh_loc,pfftIntEdt.p{1}*velocity^2)
    %legend(hca,{'\phi_E','B_{||}*scale','\int Ek dt*v'})  
    %hca.YLabel.String = 'Power [V^2/Hz]';
    %hca.XLabel.String = 'f/f_LH';
    
    hca = subplot(nrows,1,2);
    %semilogy(hca,pfftPhiE.f/flh_loc,pfftPhiE.p{1},...
    %             pfftIntEdt.f/flh_loc,pfftIntEdt.p{1}*velocity^2,...
    %             pfftBz.f/flh_loc,pfftBz.p{1}*1e-18*EBscale^2,...
    %             pfftPhiB.f/flh_loc,pfftPhiB.p{1})
    loglog(hca,pfftBz.f/flh_loc,pfftBz.p{1}*1e-18*EBscale^2,...
               pfftIntEdt.f/flh_loc,pfftIntEdt.p{1}*velocity^2)
    hca.XLim = [10^-2 2*10^1];
    legend(hca,{'B_{||}\times(B_0/n_ee\mu_0)','\int E_k dt\times v'},'location','best')
    hca.YLabel.String = '[V^2/Hz]';
    hca.XLabel.String = 'f/f_LH';
    
    hca = subplot(nrows,1,3);
    loglog(hca,wavenumber*di_loc,pfftEk.p{1},...
               wavenumber*di_loc,pfftBz.p{1},...
               wavenumber*di_loc,pfftEk.p{1}./(pfftEk.f.^2'))
    hca.XLim = klims;       
    legend(hca,{'E_k','B_{||}','E_k/f'},'location','best')%,'location','best')
    hca.YLabel.String = 'Power';
    hca.XLabel.String = 'kd_i';
    hold(hca,'on')
    loglog(xx1,yfit1,'k-',xx2,yfit2,'k-');
    %loglog(wavenumber*di_loc,yfit1fun(wavenumber*di_loc),'k-',...
    %       wavenumber*di_loc,yfit2fun(wavenumber*di_loc),'k-');
    hold(hca,'off')
    
    hca = subplot(nrows,1,4);
    loglog(hca,wavenumber*di_loc,pfftBz.p{1}*1e-18*EBscale^2,...
               pfftIntEdt.f*2*pi/velocity*di_loc,pfftIntEdt.p{1}*velocity^2,...
               di_loc/ri_loc*[1 1],hca.YLim,'k',... 
               di_loc/sqrt(ri_loc*re_loc)*[1 1],hca.YLim,'g',...
               di_loc/re_loc*[1 1],hca.YLim,'magenta')
    
    hca.XLim = klims;       
    legend(hca,{'B_{||}\times(B_0/n_ee\mu_0)','\int E_k dt\times v','k\rho_i=1','k(\rho_e\rho_i)^{1/2}=1','k\rho_e=1'},'location','best')
    hca.YLabel.String = '[V^2 /Hz]';
    hca.XLabel.String = 'kd_i';
    
end
if 0
    %%
    if 0
    irf_match_phibe_vis(hca,'best',phi_E(:,[1 1+i_v]),phi_B,v(i_v),direction,corr_dir,flim,flh_loc);
    
    irf_legend(h(5),['E_n = ' num2str(avEn,'%.f') ' mV/m   E_k = ' num2str(avEk,'%.f') ' mV/m   B_0 = ' num2str(B0,'%.f') ' nT   n = ' num2str(n_loc,'%.1f') ' cc' ],[0.02 -.5])
    irf_legend(h(5),['v_{ExB} = [' num2str(irf_norm(vExB),'%.2f') ']\times' num2str(irf_abs(vExB,1),'%.f') 'km/s' ],[0.02 -.75])
    irf_legend(h(5),['v-v_{ExB} = [' num2str(irf_norm(direction*velocity-vExB),'%.2f') ']\times' num2str(irf_abs(direction*velocity-vExB,1),'%.f') ' km/s' ],[0.02 -1.0])
    irf_legend(h(5),['\theta_{v>v_{ExB}} = ' num2str(acosd(direction*irf_norm(vExB)'),'%.f') ' ^{\circ}'],[0.02 -1.25])
    %text(h(5).XLim(1),h(5).YLim(2)*1,'a')
    delete(h(6))
    irf_plot_zoomin_lines_between_panels(h(3),h(5))
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
    end
    %%
    if doPrint; cn.print([ 'C' num2str(sc) '_' csys '_inclov_' irf_time(tint(1),'epoch>utc') '_' num2str(f_highpass,'%.1f') '_best'],'path',savePath); end
end
if 0 % Plot spectra
    %%
    flim = 0.1;
    tool.single_event  
    Ekk = irf_dot(gsmE3,direction);
    
    tint_fft = tint+5*[-1 1];
    locE = irf_tlim(irf_abs(gsmE3),tint_fft);
    locPhi = irf_abs(irf_integrate(irf_tlim(gsmE3,tint_fft)));
    locPhik = irf_abs(irf_integrate(irf_dot(irf_tlim(gsmE3,tint_fft),direction)));
    locB = irf_tlim(irf_abs(gsmB3),tint_fft);
    pfftE = irf_powerfft(locE,size(locE,1),450,0.5);
    pfftPhi = irf_powerfft(locPhi,size(locE,1),450,0.5);
    pfftPhik = irf_powerfft(locPhik,size(locPhik,1),450,0.5);
    pfftB = irf_powerfft(locB,size(locB,1),450,0.5);
    EBscale = B0.^2*1e-18/(n_loc*1e6*e*4*pi*1e-7);
    hca = subplot(2,1,1);
    semilogy(hca,pfftE.f/flh_loc,pfftE.p{4},...
                 pfftB.f/flh_loc,pfftB.p{4},...
                 pfftPhi.f/flh_loc,pfftPhi.p{4}/EBscale*velocity*1e3,...
                 pfftPhi.f/flh_loc,pfftPhi.p{4},...
                 pfftPhik.f/flh_loc,pfftPhik.p{1})
    hca = subplot(2,1,2);
    semilogy(hca,pfftE.f/flh_loc,smooth(pfftE.p{4}),...
                 pfftB.f/flh_loc,smooth(pfftB.p{4}),...
                 pfftPhi.f/flh_loc,smooth(pfftPhi.p{4}/EBscale*velocity*1e3),...
                 pfftPhi.f/flh_loc,smooth(pfftPhi.p{4}))
    %%

             %%
    hca = subplot(2,1,2);
    semilogy(hca,pfftE.f/flh_loc,smooth(pfftE.p{4}),...
                 pfftB.f/flh_loc,smooth(pfftB.p{4}),...
                 pfftPhi.f/flh_loc,smooth(pfftPhi.p{4}/EBscale*velocity*1e3),...
                 pfftPhi.f/flh_loc,smooth(pfftPhi.p{4}))
end
if 0
    %%
    ffilt = 0.01*flh_loc;
    filteroder = 3;
    Bellip = irf_filt(locBz,ffilt,0,450,filterorder);
    Bbutter = cn_filt(locBz,ffilt,13);    
    Bbutterfft = irf_powerfft(Bbutter,size(Ebutter,1),450,0.5);
    Bellipfft = irf_powerfft(Bellip,size(Ebutter,1),450,0.5);
    Bbutterfft = cn_smooth_fft(Bbutterfft,100);
    Bellipfft = cn_smooth_fft(Bellipfft,100);
    
    Eellip = irf_filt(locEk,ffilt,0,450,filterorder);
    Ebutter = cn_filt(locEk,ffilt,13);    
    Ebutterfft = irf_powerfft(Ebutter,size(Ebutter,1),450,0.5);
    Eellipfft = irf_powerfft(Eellip,size(Ebutter,1),450,0.5);
    Ebutterfft = cn_smooth_fft(Ebutterfft,100);
    Eellipfft = cn_smooth_fft(Eellipfft,100);
    
    
    loglog(Eellipfft.f,Eellipfft.p{1}*phiEscale,Ebutterfft.f,Ebutterfft.p{1}*phiEscale,...
           Bellipfft.f,Bellipfft.p{1}*1e-9*phiBscale,Bbutterfft.f,Bbutterfft.p{1}*1e-9*phiBscale)
    
    
    
    
end