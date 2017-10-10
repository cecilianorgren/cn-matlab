% Definte time interval for spectra
tint = toepoch([2008 04 22 18 07 33; 2008 04 22 18 07 36.9])';
c_eval('R_loc = irf_resamp(gsmR?,tint(1),''nearest'');',sc); R_loc = R_loc(2:4);
c_eval('Vi_loc = irf_resamp(gsmVi?,tint(1),''nearest'');',sc); Vi_loc = Vi_loc(2:4);
c_eval('Ti_loc = irf_resamp(Ti?,tint(1),''nearest'');',sc); Ti_loc = Ti_loc(2);
c_eval('Te_loc = irf_resamp(parTe?,tint(1),''nearest'');',sc); Te_loc = Te_loc(2);


    %% Some input
    flim = 0.1; % for correlation and velocity

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
    tint_fft = tint;%tint_ov_thinner2;    
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
    
    %figure(44);irf_plot({locBz,BDC})
    
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
    idxBnoise = find(plotPhiBf>1000);
    plotPhiBp(idxBnoise) = [];
    plotPhiBf(idxBnoise) = [];
    
    % remove frequencies which are below 0.03 Hz 
    lowestf = 0;0.03;
    idxBlow = find(plotPhiBf<lowestf);
    plotPhiBp(idxBlow) = [];
    plotPhiBf(idxBlow) = [];
    idxElow = find(plotPhiEf<lowestf);
    plotPhiEp(idxElow) = [];
    plotPhiEf(idxElow) = [];
    
    % plot figure
    hca = subplot(1,1,1);
    lines = loglog(hca,plotPhiBf*2*pi/velocity*re_loc,plotPhiBp,'-',...
                       plotPhiEf*2*pi/velocity*re_loc,plotPhiEp,'-');
    lines(1).Color = colors_phi(2,:);
    lines(2).Color = colors_phi(1,:);
    
    hca.YLabel.String = 'Power [V^2 /Hz]';
   %hca.YLabel.String = 'Power Spectral Density [V^2 /Hz]';
    hca.XLabel.String = 'k_{\perp}\rho_e';
    hca.Title.String = tint_str;
    %hca.XLim = [3e-2 1e3];
    %hca.YLim = [1e-6 1e10];
    tickstep = 1;
    %hca.YTick = 10.^(log10(hca.YLim(1)):tickstep:log10(hca.YLim(2)));
    hca.YTick = 10.^(-16:2:16);
    
    hold(hca,'on')
    lines = loglog(re_loc/ri_loc*[1 1],hca.YLim,... 
                   re_loc/sqrt(ri_loc*re_loc)*[1 1],hca.YLim);
    lines(1).Color = lengthcolors(2,:);
    lines(2).Color = lengthcolors(3,:);
    %lines(3).Color = lengthcolors(4,:);
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

