% test_spectra
% event20070902_loaddata.m
% run first figure of event20070902_figures.m
if 1 % See why spectra is now as it should        
    % Get propagation vector for one given highpass frequency and short     
    % tint, then use direction and velocity to integrate E again with or 
    % without filtering.
    tint_ov = toepoch([2007 09 02 15 46 30; 2007 09 02 15 48 50])';
    tint_ov1 = toepoch([2007 09 02 15 46 30; 2007 09 02 15 47 36])';
    tint_ov2 = toepoch([2007 09 02 15 47 41; 2007 09 02 15 48 50])';
    tint = toepoch([2007 09 02 15 47 45.3; 2007 09 02 15 47 46.7])';        
    sc = 3;
    csys = 'gsm';
    flim = 0.1; % *flh, for correlation and velocity
    tool.single_event  
    disp(['C = ' num2str(corr_dir(i_dir)) ',  v = ' num2str(velocity,'%.0f') ' x [' num2str(direction) '] km/s'])
    
    %% Some input
    % Detalis for filtering, plotting and comparing
    doSmooth = 1;
    v_factor = 1; % multiply v with this to adjust for other time intervals than tint
    fn = 100; % number of f intervals for smoothing  
    tint_fft = tint_ov2;%_ov; % interval to make fft on
    ffilt = 0.1*flh_loc; % for plotting spectra
    ffilt = 0.1*1/diff(tint); % make net potential difference zero
    filterorder = 3; % elliptic filter order
    
    
    
    tint_str= [datestr(irf_time(tint_fft(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint_fft(2),'epoch>datenum'),'hh:mm:ss.fff')];
    phiBscale = B0*1e-9/(n_loc*1e6*e*4*pi*1e-7);
    phiEscale = velocity*v_factor;
    
    % Prepare fields            
    % Take away NaNs so spectra can be made for longer time interval, this
    % will introduce a discontinuity, but hopefully it will not add much to
    % the power
    Ekk = irf_dot(gsmE3(~isnan(gsmE3(:,2)),:),direction); 
    maxEkk = max(Ekk(:,1));
    Bzz = irf_dot(gsmB3,z(i_dir,:));    
    maxBzz = max(Bzz(:,1));
    
    locEk = irf_tlim(Ekk,tint_fft);
    locEkFilt = irf_filt(locEk,ffilt,0,450,filterorder);
    locIntEdt = irf_integrate(locEk);  
    locIntEdtFilt = irf_integrate(locEkFilt);  
    locIntEdtFilt2 = irf_filt(locIntEdtFilt,ffilt,0,450,filterorder);  % integrated first then filtered
    
    locBz = irf_tlim(Bzz,tint_fft);
    locBzFilt = irf_filt(locBz,ffilt,0,450,filterorder);
    %locIntBdtFilt = irf_integrate(locBzFilt);  
    %locIntBdt = irf_integrate(locBz);  
        
    % Make different spectra    
    
    % Unfiltered fields
    pfftEk = irf_powerfft(locEk,size(locEk,1),450,0.5);
    pfftIntEdt = irf_powerfft(locIntEdt,size(locIntEdt,1),450,0.5);
    pfftBz = irf_powerfft(locBz,size(locBz,1),450,0.5);
    
    % Highpass filtered
    pfftIntEdtFilt = irf_powerfft(locIntEdtFilt,size(locIntEdtFilt,1),450,0.5);
    pfftBzFilt = irf_powerfft(locBzFilt,size(locBzFilt,1),450,0.5);
    %pfftIntBdtFilt = irf_powerfft(locIntBdtFilt,size(locBzFilt,1),450,0.5);
    %pfftIntBdt = irf_powerfft(locIntBdt,size(locBz,1),450,0.5);
    
    if doSmooth % make spectra smooth
        %%
        disp(['Smoothing spectra. Using ' num2str(fn) ' frequency intervals.'])
        %pfftIntBdtFilt = cn_smooth_fft(pfftIntBdtFilt,fn);
        %pfftIntBdt = cn_smooth_fft(pfftIntBdt,fn);
        pfftIntEdtFilt = cn_smooth_fft(pfftIntEdtFilt,fn);
        pfftIntEdt = cn_smooth_fft(pfftIntEdt,fn);
        pfftEk = cn_smooth_fft(pfftEk,fn);
        pfftBz = cn_smooth_fft(pfftBz,fn);
        pfftBzFilt = cn_smooth_fft(pfftBzFilt,fn);
    end    
    
    if 0 % plot wavelet
        %%
        units =irf_units;
        wavewidth = 1*5.36;
        figure(76)
        waveEk = irf_wavelet(locEk,'wavelet_width',wavewidth);
        waveIntEdt = irf_wavelet(locIntEdt,'wavelet_width',wavewidth);
        waveBz = irf_wavelet(locBz,'wavelet_width',wavewidth);
        waveEB = waveEk;
        waveEB.p{1} = sqrt(waveEk.p{1}*1e-6./(waveBz.p{1}*1e-18));
        
        hca = subplot(3,1,1);
        pcolor(hca,waveIntEdt.t,waveIntEdt.f,log10(waveIntEdt.p{1}'*phiEscale^2))        
        shading flat;
        hca.YScale = 'log';
        hca.YLabel.String = '\phi_E^2 [V^2/Hz]';
        
        hca = subplot(3,1,2);
        pcolor(hca,waveBz.t,waveBz.f,log10(waveBz.p{1}'*phiBscale^2))        
        shading flat;
        hca.YScale = 'log';
        hca.YLabel.String = '\phi_B^2 [V^2/Hz]';
        
        hca = subplot(3,1,3);
        pcolor(hca,waveEB.t,waveEB.f,log10(waveEB.p{1}'/units.c))        
        shading flat;
        hca.YScale = 'log';
        hca.YLabel.String = 'log_{10} E/B/c';
        hc = colorbar;
        hca.CLim = [-2 0];
    end       
    if 0
        loglog(pfftIntBdtFilt.f,pfftIntBdtFilt.p{1},pfftIntBdt.f,pfftIntBdt.p{1},...
           pfftIntEdtFilt.f,pfftIntEdtFilt.p{1},pfftIntEdt.f,pfftIntEdt.p{1},...
           pfftEkFilt.f,pfftEkFilt.p{1},pfftIntEdt.f,pfftIntEdt.p{1});
        legend('Bfilt','B','Efilt','E')
    end
    
    
    % Plot different spectra
    fig = figure(98);
    
    xticks = 10.^[-1 0 1 2];
    xlim = [1/diff(tint_fft) 250];
    
    nPanels = 10;
    for ii = 1:nPanels; h(ii) = subplot(nPanels/2,2,ii); end    
    isub = 1;
    if 1 % Ek
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftEk.f,pfftEk.p{1});
        hca.XLim = xlim;
        hca.XTick = xticks;
        hca.YLabel.String = 'Power [(mV/m)^2/Hz]';
        hca.Title.String = 'unfiltered';
        
        hca = h(isub); isub = isub + 1;
        irf_plot(hca,locEk);
        hca.YLabel.String = 'E [mV/m]';
        hca.Title.String = 'unfiltered';
        irf_zoom(hca,'x',tint_fft);        
    end
    if 1 % Bz
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftBz.f,pfftBz.p{1}); 
        hca.XTick = xticks;
        hca.XLim = xlim;
        hca.YLabel.String = 'Power [(nT)^2/Hz]';
        hca.Title.String = 'unfiltered';
        
        hca = h(isub); isub = isub + 1;
        irf_plot(hca,locBz);
        hca.YLabel.String = 'B [nT]';
        hca.Title.String = 'unfiltered';
        irf_zoom(hca,'x',tint_fft);
    end
    
    if 1 % int E dt
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftIntEdt.f,pfftIntEdt.p{1});
        hca.XLim = xlim;
        hca.XTick = 10.^[-1 0 1 2];
        hca.YLabel.String = 'Power [(mV/m*s)^2/Hz]';
        hca.Title.String = 'unfiltered';
        
        hca = h(isub); isub = isub + 1;
        irf_plot(hca,locIntEdt);
        hca.Title.String = 'unfiltered';
        hca.YLabel.String = '\int E dt [mV/m*s]';
        irf_zoom(hca,'x',tint_fft);
    end
    
    if 1 % filtered int E dt
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftIntEdtFilt.f,pfftIntEdtFilt.p{1});
        hca.XLim = xlim;
        hca.XTick = 10.^[-1 0 1 2];
        hca.YLabel.String = 'Power [(mV/m*s)^2/Hz]';
        hca.Title.String = ['highpass filtered at f=' num2str(ffilt,'%.2f') ' Hz'];
        
        hca = h(isub); isub = isub + 1;
        irf_plot(hca,locIntEdtFilt);
        hca.Title.String = ['highpass filtered at f=' num2str(ffilt,'%.2f') ' Hz'];
        irf_zoom(hca,'x',tint_fft);
    end
    if 1 % compare filtered int E dt, unfiltered int E dt and Bz, all in units of V    
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftIntEdtFilt.f,pfftIntEdtFilt.p{1}*phiEscale^2,...
                   pfftIntEdt.f,pfftIntEdt.p{1}*phiEscale^2,...
                   pfftBzFilt.f,pfftBzFilt.p{1}*1e-18*phiBscale^2,...
                   pfftBz.f,pfftBz.p{1}*1e-18*phiBscale^2);
        hca.YLabel.String = 'Power [V^2/Hz]';
        irf_legend(hca,{'\phi_{E,filt}','\phi_{E}','\phi_{B,filt}','\phi_{B}'},[0.98 0.95])
        hca.XLim = xlim;
        hca.YLim = 10.^[-10 10];
        hca.YTick = 10.^[-10:5:10];
        hca.XTick = xticks;
        %hca.YLim = 10.^[-5 5];
    end
    if 1 % compare filtered int E dt, unfiltered int E dt and Bz, all in units of V    
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftIntEdtFilt.f*2*pi*di_loc/velocity,pfftIntEdtFilt.p{1}*phiEscale^2,...
                   pfftIntEdt.f*2*pi*di_loc/velocity,pfftIntEdt.p{1}*phiEscale^2,...
                   pfftBzFilt.f*2*pi*di_loc/velocity,pfftBzFilt.p{1}*1e-18*phiBscale^2,...
                   pfftBz.f*2*pi*di_loc/velocity,pfftBz.p{1}*1e-18*phiBscale^2);
        hca.YLabel.String = 'Power [V^2/Hz]';
        irf_legend(hca,{'\phi_{E,filt}','\phi_{E}','\phi_{B,filt}','\phi_{B}'},[0.98 0.95])
        hca.XLim = xlim*2*pi*di_loc/velocity;
        hca.XTick = 10.^[-1 0 1 2 3 4];
        hca.YLim = 10.^[-10 10];
        hca.YTick = 10.^[-10:5:10];
        hca.XLabel.String = 'kd_i';
    end
    %%      
    if 1 % compare filtered int E dt, unfiltered int E dt and Bz, all in units of V    
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftIntEdtFilt.f,pfftIntEdtFilt.p{1}*phiEscale^2,...
                   pfftIntEdt.f,pfftIntEdt.p{1}*phiEscale^2,...
                   pfftBzFilt.f,pfftBzFilt.p{1}*1e-18*phiBscale^2,...
                   pfftBz.f,pfftBz.p{1}*1e-18*phiBscale^2);
        legend(hca,'\phi_{E,filt}','\phi_{E}','\phi_{B,filt}','\phi_{B}','location','best')
        hca.XLim = [6e-1 250];                
    end
    %%
    %lambda = velocity./pfftEk.f; % km    
    %wavenumber = 2*pi./lambda; % 1/km
    %wavenumber = pfftEk.f*2*pi/velocity;

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
if 1 % See why spectra is now as it should        
    % Get propagation vector for one given highpass frequency and short     
    % tint, then use direction and velocity to integrate E again with or 
    % without filtering.
    tint_ov = toepoch([2007 09 02 15 46 30; 2007 09 02 15 48 50])';
    tint_ov1 = toepoch([2007 09 02 15 46 30; 2007 09 02 15 47 36])';
    tint_ov2 = toepoch([2007 09 02 15 47 41; 2007 09 02 15 48 50])';
    tint = toepoch([2007 09 02 15 47 45.3; 2007 09 02 15 47 46.7])';        
    sc = 3;
    csys = 'gsm';
    flim = 0.1; % *flh, for correlation and velocity
    tool.single_event  
    disp(['C = ' num2str(corr_dir(i_dir)) ',  v = ' num2str(velocity,'%.0f') ' x [' num2str(direction) '] km/s'])
    
    %% Some input
    % Detalis for filtering, plotting and comparing
    doSmooth = 1;
    v_factor = 1; % multiply v with this to adjust for other time intervals than tint
    fn = 100; % number of f intervals for smoothing  
    tint_fft = tint_ov;%_ov; % interval to make fft on
    ffilt = 0.1*flh_loc; % for plotting spectra
    ffilt = 0.1*1/diff(tint); % make net potential difference zero
    filterorder = 3; % elliptic filter order
    
    
    
    tint_str= [datestr(irf_time(tint_fft(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint_fft(2),'epoch>datenum'),'hh:mm:ss.fff')];
    phiBscale = B0*1e-9/(n_loc*1e6*e*4*pi*1e-7);
    phiEscale = velocity*v_factor;
    
    % Prepare fields            
    % Take away NaNs so spectra can be made for longer time interval, this
    % will introduce a discontinuity, but hopefully it will not add much to
    % the power
    %Ekk = irf_dot(gsmE3(~isnan(gsmE3(:,2)),:),direction); 
    %maxEkk = max(Ekk(:,1));
    Bzz = irf_dot(gsmB3,z(i_dir,:));    
    maxBzz = max(Bzz(:,1));
    locE = irf_tlim(gsmE3,tint_fft);
    
    locEkNanIncl = locE(:,[1 2]); locE(:,2) = locE(:,2:4)*direction';
    EisNaN = isnan(locEkNanIncl(:,2));
    nanstart = find(EisNaN>0,1,'first');
    nanstop = find(EisNaN>0,1,'last');
    nanall = find(EisNaN>0);
    E1 = locEkNanIncl(nanstart-1,2);
    E2 = locEkNanIncl(nanstop+1,2);
    % linear interpolation
    Einterp = E1+(nanall-nanall(1))/(nanall(end)-nanall(1))*(E2-E1);
    locEkInterpNan = locEkNanIncl; locEkInterpNan(nanall,2) = Einterp;    
    
    locEk = irf_tlim(Ekk,tint_fft);
    locEkFilt = irf_filt(locEk,ffilt,0,450,filterorder);
    locIntEdt = irf_integrate(locEk);  
    locIntEdtFilt = irf_integrate(locEkFilt);  
    locIntEdtFilt2 = irf_filt(locIntEdtFilt,ffilt,0,450,filterorder);  % integrated first then filtered
    
    locBz = irf_tlim(Bzz,tint_fft);
    locBzFilt = irf_filt(locBz,ffilt,0,450,filterorder);
    %locIntBdtFilt = irf_integrate(locBzFilt);  
    %locIntBdt = irf_integrate(locBz);  
        
    % Make different spectra    
    
    % Unfiltered fields
    pfftEk = irf_powerfft(locEk,size(locEk,1),450,0.5);
    pfftIntEdt = irf_powerfft(locIntEdt,size(locIntEdt,1),450,0.5);
    pfftBz = irf_powerfft(locBz,size(locBz,1),450,0.5);
    pfftEkInterpNan = irf_powerfft(locEkInterpNan,size(locEkInterpNan,1),450,0.5);
    
    % Highpass filtered
    pfftIntEdtFilt = irf_powerfft(locIntEdtFilt,size(locIntEdtFilt,1),450,0.5);
    pfftBzFilt = irf_powerfft(locBzFilt,size(locBzFilt,1),450,0.5);
    %pfftIntBdtFilt = irf_powerfft(locIntBdtFilt,size(locBzFilt,1),450,0.5);
    %pfftIntBdt = irf_powerfft(locIntBdt,size(locBz,1),450,0.5);
    
    %loglog(pfftEk.f,pfftEk.p{1});
    %%
    if doSmooth % make spectra smooth
        %%
        disp(['Smoothing spectra. Using ' num2str(fn) ' frequency intervals.'])
        %pfftIntBdtFilt = cn_smooth_fft(pfftIntBdtFilt,fn);
        %pfftIntBdt = cn_smooth_fft(pfftIntBdt,fn);
        pfftIntEdtFilt = cn_smooth_fft(pfftIntEdtFilt,fn);
        pfftIntEdt = cn_smooth_fft(pfftIntEdt,fn);
        pfftEk = cn_smooth_fft(pfftEk,fn);
        pfftEkInterpNan = cn_smooth_fft(pfftEkInterpNan,fn);
        pfftBz = cn_smooth_fft(pfftBz,fn);
        pfftBzFilt = cn_smooth_fft(pfftBzFilt,fn);
    end    
    
    
    % Plot different spectra
    fig = figure(98);
    
    xticks = 10.^[-1 0 1 2];
    xlim = [1/diff(tint_fft) 250];
    
    nPanels = 10;
    for ii = 1:nPanels; h(ii) = subplot(nPanels/2,2,ii); end    
    isub = 1;
    if 1 % Ek
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftEk.f,pfftEk.p{1});
        hca.XLim = xlim;
        hca.XTick = xticks;
        hca.YLabel.String = 'Power [(mV/m)^2/Hz]';
        hca.Title.String = 'unfiltered';
        
        hca = h(isub); isub = isub + 1;
        irf_plot(hca,locEk);
        hca.YLabel.String = 'E [mV/m]';
        hca.Title.String = 'unfiltered';
        irf_zoom(hca,'x',tint_fft);        
    end
    if 1 % Bz
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftBz.f,pfftBz.p{1}); 
        hca.XTick = xticks;
        hca.XLim = xlim;
        hca.YLabel.String = 'Power [(nT)^2/Hz]';
        hca.Title.String = 'unfiltered';
        
        hca = h(isub); isub = isub + 1;
        irf_plot(hca,locBz);
        hca.YLabel.String = 'B [nT]';
        hca.Title.String = 'unfiltered';
        irf_zoom(hca,'x',tint_fft);
    end
    
    if 1 % int E dt
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftIntEdt.f,pfftIntEdt.p{1});
        hca.XLim = xlim;
        hca.XTick = 10.^[-1 0 1 2];
        hca.YLabel.String = 'Power [(mV/m*s)^2/Hz]';
        hca.Title.String = 'unfiltered';
        
        hca = h(isub); isub = isub + 1;
        irf_plot(hca,locIntEdt);
        hca.Title.String = 'unfiltered';
        hca.YLabel.String = '\int E dt [mV/m*s]';
        irf_zoom(hca,'x',tint_fft);
    end
    
    if 1 % filtered int E dt
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftIntEdtFilt.f,pfftIntEdtFilt.p{1});
        hca.XLim = xlim;
        hca.XTick = 10.^[-1 0 1 2];
        hca.YLabel.String = 'Power [(mV/m*s)^2/Hz]';
        hca.Title.String = ['highpass filtered at f=' num2str(ffilt,'%.2f') ' Hz'];
        
        hca = h(isub); isub = isub + 1;
        irf_plot(hca,locIntEdtFilt);
        hca.Title.String = ['highpass filtered at f=' num2str(ffilt,'%.2f') ' Hz'];
        irf_zoom(hca,'x',tint_fft);
    end
    if 1 % compare filtered int E dt, unfiltered int E dt and Bz, all in units of V    
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftIntEdtFilt.f,pfftIntEdtFilt.p{1}*phiEscale^2,...
                   pfftEkInterpNan.f,pfftEkInterpNan.p{1}*phiEscale^2,...
                   pfftIntEdt.f,pfftIntEdt.p{1}*phiEscale^2,...
                   pfftBzFilt.f,pfftBzFilt.p{1}*1e-18*phiBscale^2,...
                   pfftBz.f,pfftBz.p{1}*1e-18*phiBscale^2);
        hca.YLabel.String = 'Power [V^2/Hz]';
        irf_legend(hca,{'\phi_{E,filt}','\phi_{E,interpNaNs}','\phi_{E}','\phi_{B,filt}','\phi_{B}'},[0.98 0.95])
        hca.XLim = xlim;
        hca.YLim = 10.^[-10 10];
        hca.YTick = 10.^[-10:5:10];
        hca.XTick = xticks;
        %hca.YLim = 10.^[-5 5];
    end
    if 1 % compare filtered int E dt, unfiltered int E dt and Bz, all in units of V    
        hca = h(isub); isub = isub + 1;
        loglog(hca,pfftIntEdtFilt.f*2*pi*di_loc/velocity,pfftIntEdtFilt.p{1}*phiEscale^2,...
                   pfftIntEdt.f*2*pi*di_loc/velocity,pfftIntEdt.p{1}*phiEscale^2,...
                   pfftBzFilt.f*2*pi*di_loc/velocity,pfftBzFilt.p{1}*1e-18*phiBscale^2,...
                   pfftBz.f*2*pi*di_loc/velocity,pfftBz.p{1}*1e-18*phiBscale^2);
        hca.YLabel.String = 'Power [V^2/Hz]';
        irf_legend(hca,{'\phi_{E,filt}','\phi_{E}','\phi_{B,filt}','\phi_{B}'},[0.98 0.95])
        hca.XLim = xlim*2*pi*di_loc/velocity;
        hca.XTick = 10.^[-1 0 1 2 3 4];
        hca.YLim = 10.^[-10 10];
        hca.YTick = 10.^[-10:5:10];
        hca.XLabel.String = 'kd_i';
    end
    if 0
        %%
        pfftEkInterpNan = cn_smooth_fft(pfftEkInterpNan,fn);
        loglog(pfftEk.f,pfftEk.p{1},pfftEkInterpNan.f,pfftEkInterpNan.p{1});
        legend('cut out nans','interpolated nans')
    end
end
if 1 % different nan treatments
    %%
    % Get propagation vector for one given highpass frequency and short     
    % tint, then use direction and velocity to integrate E again with or 
    % without filtering.
    tint_ov = toepoch([2007 09 02 15 46 30; 2007 09 02 15 48 50])';
    tint_ov1 = toepoch([2007 09 02 15 46 30; 2007 09 02 15 47 36])';
    tint_ov2 = toepoch([2007 09 02 15 47 41; 2007 09 02 15 48 50])';
    tint = toepoch([2007 09 02 15 47 45.3; 2007 09 02 15 47 46.7])';        
    sc = 3;
    csys = 'gsm';
    flim = 0.1; % *flh, for correlation and velocity
    tool.single_event  
    disp(['C = ' num2str(corr_dir(i_dir)) ',  v = ' num2str(velocity,'%.0f') ' x [' num2str(direction) '] km/s'])
    
    %% Some input
    % Detalis for filtering, plotting and comparing
    doSmooth = 1;
    v_factor = 1; % multiply v with this to adjust for other time intervals than tint
    fn = 100; % number of f intervals for smoothing  
    tint_fft = tint_ov;%_ov; % interval to make fft on
    ffilt = 0.1*flh_loc; % for plotting spectra
    ffilt = 0.1*1/diff(tint); % make net potential difference zero
    filterorder = 3; % elliptic filter order
    
    
    
    tint_str= [datestr(irf_time(tint_fft(1),'epoch>datenum'),'hh:mm:ss.fff') ' - ' datestr(irf_time(tint_fft(2),'epoch>datenum'),'hh:mm:ss.fff')];
    phiBscale = B0*1e-9/(n_loc*1e6*e*4*pi*1e-7);
    phiEscale = velocity*v_factor;
    
    % Prepare fields            
    % Take away NaNs so spectra can be made for longer time interval, this
    % will introduce a discontinuity, but hopefully it will not add much to
    % the power
    
    Bzz = irf_dot(gsmB3,z(i_dir,:));    
    maxBzz = max(Bzz(:,1));
    locE = irf_tlim(gsmE3,tint_fft);
    locEk = locE(:,1:2); locEk(:,2) = locE(:,2:4)*direction'; % includes NaNs
    maxEk = max(locEk(:,1));
    
    % find nans
    EisNaN = isnan(locEk(:,2));
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
    locIntEkNanInterp = irf_integrate(locEkNanInterp);  
    pfftIntEkNanInterp = irf_powerfft(locIntEkNanInterp,size(locIntEkNanInterp,1),450,0.5);
    locIntEkNanInterpFilt = irf_filt(locIntEkNanInterp,ffilt,0,450,filterorder);
    pfftIntEkNanInterpFilt = irf_powerfft(locIntEkNanInterpFilt,size(locIntEkNanInterpFilt,1),450,0.5);
    
    % set all nans to zero
    locEkNanZero = locEk(:,[1 2]);
    locEkNanZero(nanall,2) = locEkNanZero(nanall,2)*0;
    pfftEkNanZero = irf_powerfft(locEkNanZero,size(locEkNanZero,1),450,0.5);
    locIntEkNanZero = irf_integrate(locEkNanZero);  
    pfftIntEkNanZero = irf_powerfft(locIntEkNanZero,size(locIntEkNanZero,1),450,0.5);    
    locIntEkNanZeroFilt = irf_filt(locIntEkNanZero,ffilt,0,450,filterorder);
    pfftIntEkNanZeroFilt = irf_powerfft(locIntEkNanZeroFilt,size(locIntEkNanZeroFilt,1),450,0.5);
    
    % cut off from time series
    locEkNanCut = locEk; locEkNanCut(nanall,:) = [];
    pfftEkNanCut = irf_powerfft(locEkNanCut,size(locEkNanCut,1),450,0.5);
    locIntEkNanCut = irf_integrate(locEkNanCut);  
    pfftIntEkNanCut = irf_powerfft(locIntEkNanCut,size(locIntEkNanCut,1),450,0.5);
    locIntEkNanCutFilt = irf_filt(locIntEkNanCut,ffilt,0,450,filterorder);
    pfftIntEkNanCutFilt = irf_powerfft(locIntEkNanCutFilt,size(locIntEkNanCutFilt,1),450,0.5);
    %if doSmooth, pfftIntEkNanCut = cn_smooth_fft(pfftIntEkNanCut); pfftEkNanCut = cn_smooth_fft(pfftEkNanCut); end
    
    % highpass filter fields
    locEkFilt = irf_filt(locEk,ffilt,0,450,filterorder);
    
    locIntEdtFilt = irf_integrate(locEkFilt);  
    locFiltIntEdt = irf_filt(locIntEdtFilt,ffilt,0,450,filterorder);  % integrated first then filtered
    
    locBz = irf_tlim(Bzz,tint_fft);
    locBzFilt = irf_filt(locBz,ffilt,0,450,filterorder);
    
    % Unfiltered fields
    locIntEdt = irf_integrate(locEk);  
    
    pfftEk = irf_powerfft(locEk,size(locEk,1),450,0.5);
    pfftIntEdt = irf_powerfft(locIntEdt,size(locIntEdt,1),450,0.5);
    pfftBz = irf_powerfft(locBz,size(locBz,1),450,0.5);    
    
    % Highpass filtered
    pfftIntEdtFilt = irf_powerfft(locIntEdtFilt,size(locIntEdtFilt,1),450,0.5);
    pfftBzFilt = irf_powerfft(locBzFilt,size(locBzFilt,1),450,0.5);

    % Smoothen spectra
    if doSmooth
        clear pffts        
        pffts = who('pfft*') 
        for ii = 1:numel(pffts)
            disp(['Smoothing ' pffts{ii}])
            eval([pffts{ii} ' = cn_smooth_fft(' pffts{ii} ',fn);'])
        end
    end
    
    if 0;doSmooth % make spectra smooth
        %%
        disp(['Smoothing spectra. Using ' num2str(fn) ' frequency intervals.'])
        %pfftIntBdtFilt = cn_smooth_fft(pfftIntBdtFilt,fn);
        %pfftIntBdt = cn_smooth_fft(pfftIntBdt,fn);
        pfftIntEdtFilt = cn_smooth_fft(pfftIntEdtFilt,fn);
        pfftIntEdt = cn_smooth_fft(pfftIntEdt,fn);
        pfftEk = cn_smooth_fft(pfftEk,fn);
        pfftEkInterpNan = cn_smooth_fft(pfftEkInterpNan,fn);
        pfftBz = cn_smooth_fft(pfftBz,fn);
        pfftBzFilt = cn_smooth_fft(pfftBzFilt,fn);
    end    
    if 1 % compare filtered int E dt, unfiltered int E dt and Bz, all in units of V    
        
        hca = subplot(1,1,1);hold(hca,'on');
        colors = hca.ColorOrder;        
        plot(hca,pfftIntEkNanCut.f,pfftIntEkNanCut.p{1}*phiEscale^2,'color',colors(1,:))
        plot(hca,pfftIntEkNanCutFilt.f,pfftIntEkNanCutFilt.p{1}*phiEscale^2,'color',colors(1,:),'linestyle','--')
        plot(hca,pfftIntEkNanInterp.f,pfftIntEkNanInterp.p{1}*phiEscale^2,'color',colors(2,:))
        plot(hca,pfftIntEkNanInterpFilt.f,pfftIntEkNanInterpFilt.p{1}*phiEscale^2,'color',colors(2,:),'linestyle','--')
        plot(hca,pfftIntEkNanZero.f,pfftIntEkNanZero.p{1}*phiEscale^2,'color',colors(3,:))                                     
        plot(hca,pfftIntEkNanZeroFilt.f,pfftIntEkNanZeroFilt.p{1}*phiEscale^2,'color',colors(3,:),'linestyle','--')
        plot(hca,pfftBz.f,pfftBz.p{1}*1e-18*phiBscale^2,'color',[0 0 0])
        plot(hca,pfftBzFilt.f,pfftBzFilt.p{1}*1e-18*phiBscale^2,'color',[0 0 0],'linestyle','--');     
        hold(hca,'off');
        hca.Title.String = 'Effect on power spectra for different treatments of NaNs';
        hca.YLabel.String = 'Power [V^2/Hz]';
        hca.XLabel.String = 'f';
        %irf_legend(hca,{'\phi_{E,nancut}','\phi_{E,naninterp}','\phi_{E,nanzero}','\phi_{E,nancut,filt}','\phi_{E,naninterp,filt}','\phi_{E,nanzero,filt}','\phi_{B,filt}','\phi_{B}'},[0.98 0.95])
        legend(hca,{'\phi_{E,nancut}','\phi_{E,nancut,filt}','\phi_{E,naninterp}','\phi_{E,naninterp,filt}','\phi_{E,nanzero}','\phi_{E,nanzero,filt}','\phi_{B}','\phi_{B,filt}'})
        hca.XLim = [2e-2 300];
        hca.YLim = 10.^[-10 10];
        hca.YTick = 10.^(-10:5:10);
        hca.XTick = 10.^(-20:20);
        hca.YLim = 10.^[-5 9];
        hca.YScale = 'log';
        hca.XScale = 'log';
        hca.Box = 'on';
    end
    if 1
        figure(65)
    end
end
if 0 % compare power spectrograms of discrete and less discrete signals
    %%
    
    tsGrow = locEk; tsGrow(:,2) = (tsGrow(:,1)-tsGrow(1,1))/(tsGrow(end,1)-tsGrow(1,1))*30;
    pfftGrow = irf_powerfft(tsGrow,size(tsGrow,1),450,0.5);
    
    nT = size(tsDisc,1); nT1 = fix(nT/2);
    slopewidth = 0.05; % percentage
    startT = fix(nT/2-nT*slopewidth/2); stopT = fix(nT/2+nT*slopewidth/2);
    tsFlatGrow = locEk; 
    hmin = 5; hmax = 20;
    tsFlatGrow(1:startT,2) = hmin*ones(startT,1); tsFlatGrow((stopT+1):end,2) = hmax*ones(nT-stopT,1);
    tsFlatGrow(startT:stopT,2) = ((startT:stopT)-startT)/(stopT-startT)*(hmax-hmin)+hmin;
    pfftFlatGrow = irf_powerfft(tsFlatGrow,size(tsFlatGrow,1),450,0.5);
    
    
    tsDisc = locEk; nT = size(tsDisc,1); nT1 = fix(nT/2);
    tsDisc(1:nT1,2) = 5*ones(nT1,1); tsDisc((nT1+1):end,2) = 20*ones(nT-nT1,1);
    pfftDisc = irf_powerfft(tsDisc,size(tsDisc,1),450,0.5);
    
    tsGrowDisc = locEk; tsGrowDisc(:,2) = round(tsGrowDisc(:,1)*0.1-tsGrowDisc(1,1)*0.1)/(tsGrowDisc(end,1)-tsGrowDisc(1,1))*30*10;
    pfftGrowDisc = irf_powerfft(tsGrowDisc,size(tsGrowDisc,1),450,0.5);
    
    
    hca = subplot(2,1,1);
    irf_plot(hca,{tsDisc,tsGrowDisc,tsGrow,tsFlatGrow},'comp');
    hca.YLabel.String = 'Signal';
    hca.Title.String = 'Power spectra of different fields';
    legend('Discontinuity','Stair discontinuities','Growing field','Growing + edges','location','northeastoutside')
    
    hca = subplot(2,1,2);
    allp = [pfftDisc.p{1}',pfftGrowDisc.p{1}',pfftGrow.p{1}',pfftFlatGrow.p{1}'];
    lims = [min(min(allp)) max(max(allp))];
    loglog(hca,pfftDisc.f,allp);
    hca.YLabel.String = 'Power';
    hca.YLim = lims.*[-1 1]*10;
    
    legend('Discontinuity','Stair discontinuities','Growing field','Growing + edges','location','northeastoutside')
end