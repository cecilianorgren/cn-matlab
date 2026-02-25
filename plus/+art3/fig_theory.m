cd /Users/Cecilia/Research/EH2/Runs-solver_Buneman_kvR/
doLoad = 1;
doClearup = 1;
doPlot = 1;
doPhi = 0;
if doLoad
    loadPath = '/Users/Cecilia/Research/EH2/BeamSolver/';
    matName = '2015-03-21T012941_beam_solver_67kx27Sx20Rx4Te2x3Ti.mat';
    matName = '2015-03-25T141546_beam_solver_67kx27Sx20Rx4Te2x3Ti.mat';
    matName = '2015-03-26T122805_beam_solver_67kx27Sx20Rx4Te2x3Ti.mat';
    matName = '2015-03-27T111315_beam_solver_67kx27Sx20Rx4Te2x3Ti.mat';
    matName = '2015-03-27T192411_beam_solver_67kx27Sx20Rx2Te2x2Ti';
    matName = '2015-05-13T143716_beam_solver_67kx29Sx2Rx1Te2x1Ti.mat';
    matName = '2015-05-13T183528_beam_solver_67kx29Sx20Rx1Te2x1Ti.mat';
    matNAme = '2015-05-17T214209_beam_solver_110kx29Sx20Rx1Te2x1Ti.mat';
    
    load([loadPath matName])
    if 0
    B = 25;
    n = 0.06;
    no = 0;

    % Ions
    omega_pi = irf_plasma_calc(B,n,no,Te1,Ti,'Fpp')*2*pi; % rad/s
    vthi = irf_plasma_calc(B,n,no,Te1,Ti,'Vtp'); % m/s 
    vdi = 0;
    % Total electron plasma frequency
    omega_pe = irf_plasma_calc(B,n,no,Te1,Ti,'Fpe')*2*pi; % rad/s
    % n1= n*(1-R); % background
    % n2= n*R;     % beam 
    % Background electrons
    % omega_pe1 = irf_plasma_calc(B,n1,no,Te1,Ti,'Fpe')*2*pi; % rad/s
    vthe1 = irf_plasma_calc(B,n,no,Te1,Ti,'Vte'); % m/s 
    vde1 = 0; % m/s
    % Beam electrons
    % omega_pe2 = irf_plasma_calc(B,n2,no,Te2,Ti,'Fpe')*2*pi; % rad/s
    vthe2 = irf_plasma_calc(B,n,no,Te2,Ti,'Vte'); % m/s 
    vde2 = S*vthe1;
    % Other
    lamD = irf_plasma_calc(B,n,no,Te1,Ti,'Ld'); % m
    omega_bune = omega_pe^(1/3)*omega_pi^(2/3)*16^(-1/3);

    normal_frequency = omega_pe; % total electron plasma frequency
    normal_length = lamD; % Vte./Wpe/sqrt(2); % total density, background temperature
    normal_frequency_string = '\omega_{pe}';
    normal_length_string = '\lambda_D';
    end
end
if doClearup    
    %% clear up the matrix from phoney values
    neggrowth = find(wimax<0);
    wimax(neggrowth) = 0;
    wrmax(neggrowth) = 0;
    vphmax(neggrowth) = 0;
    
    if 0
        negwrowth = x_imag_store(x_real_store<0);
        k_store = repmat(k',1,nv,nR);
        vph_store = x_real_store/normal_frequency./k_store/sqrt(2);
        x_imag_store(x_real_store<0) = NaN;    
        vph_store(x_real_store<0) = NaN;
        x_real_store(x_real_store<0) = NaN;
        x_real_store(x_real_store/normal_frequency>2.4) = NaN;
    end

end
if doPlot               
    if 0
    vph_plot=vph_max_store;
    vph_plot(vph_max_store>2)=NaN;
    vph_plot(vph_max_store<1e-3)=NaN;
    x_imag_plot = x_imag_max_store_2Dvalues;
    x_imag_plot(x_imag_max_store_2Dvalues/normal_frequency>0.2)=NaN;
    x_imag_plot(vph_max_store>2)=NaN;
    x_imag_plot(x_imag_max_store_2Dvalues/normal_frequency<-0.02)=NaN;
    end
    if doPhi
        %%
        units = irf_units;
        vthebeam = cn_eV2v(Te2,'eV'); % km/s
        vthebg = cn_eV2v(Te1,'eV'); % km/s
        vbeam = repmat(S',1,numel(R))*vthebg;
        ephi = (vbeam-1*vthebeam-vph_plot*vthebg).^2*1e6*units.me/2/units.e/Te1;
    end
    %% make figure        
    set(0,'defaultAxesFontSize',16);
    set(0,'DefaultTextFontSize',16);
    set(0,'defaultAxesFontUnits','pixels');
    set(0,'defaultTextFontUnits','pixels');
    nPlots = 6;
    nrows = 3;
    for kk = 1:nPlots; h(kk)=subplot(nrows,2,kk); end
    isub=1;    
    klim = [0 2];
    
    % indices of what to plot
    iRplot = 5;
    iTe2plot = 1;
    iTiplot = 1;
    wilim = 8;
    wrlim = 50;
    if 1 % wi vs k,S
        hca = h(isub); isub=isub+1;           
        wiplot = squeeze(wi(:,:,iRplot,iTe2plot,iTiplot)')*sqrt(1836); % in units of ion plasma frequency
        pcolor(hca,[0 k(2:end)],S,wiplot)
        
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'S')        
        ch = colorbar('peer',hca);        
        ylabel(ch,'\omega_{i}/\omega_{pi}')  
        
        set(hca,'clim',wilim*[-1 1],'ylim',S([1 end]),'xlim',klim)
        set(ch,'ylim',wilim*[-1 1])
        shading(hca,'flat')
        %title(hca,['Growthrate at R = ' num2str(R(iRplot))])
        xlab{isub-1} = [' R = ' num2str(R(iRplot))];
        
    end
    if 1 % wr vs k,S
        hca = h(isub); isub=isub+1;           
        wrplot = squeeze(wr(:,:,iRplot,iTe2plot,iTiplot)')*sqrt(1836); % in units of ion plasma frequency
        pcolor(hca,[0 k(2:end)],S,wrplot)
        
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'S')        
        ch = colorbar('peer',hca);        
        ylabel(ch,'\omega_{r}/\omega_{pi}')  
        
        set(hca,'clim',wrlim*[-1 1],'ylim',S([1 end]),'xlim',klim)
        set(ch,'ylim',wrlim*[0 1])
        shading(hca,'flat')
        
        %title(hca,['Frequency at R = ' num2str(R(iRplot))])
        xlab{isub-1} = [' R = ' num2str(R(iRplot))];
    end
    if 0 % wi vs k,S
        hca = h(isub); isub=isub+1;           
        wiguessplot = squeeze(wiin(:,:,iRplot,iTe2plot,iTiplot)')*sqrt(1836); % in units of ion plasma frequency
        pcolor(hca,[0 k(2:end)],S,wiguessplot)
        
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'S')        
        ch = colorbar('peer',hca);        
        ylabel(ch,'\omega_{i}/\omega_{pi}')  
        
        set(hca,'clim',clim*[-1 1],'ylim',S([1 end]),'xlim',klim)
        set(ch,'ylim',clim*[-1 1])
        shading(hca,'flat')
        %title(hca,['Growthrate at R = ' num2str(R(iRplot))])
        xlab{isub-1} = [' R = ' num2str(R(iRplot))];
        
    end
    if 0 % wr vs k,S
        hca = h(isub); isub=isub+1;           
        wrguessplot = squeeze(wrin(:,:,iRplot,iTe2plot,iTiplot)')*sqrt(1836); % in units of ion plasma frequency
        pcolor(hca,[0 k(2:end)],S,wrguessplot)
        
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'S')        
        ch = colorbar('peer',hca);        
        ylabel(ch,'\omega_{r}/\omega_{pi}')  
        
        set(hca,'clim',40*[-1 1],'ylim',S([1 end]),'xlim',klim)
        set(ch,'ylim',[0 40])
        shading(hca,'flat')
        
        %title(hca,['Frequency at R = ' num2str(R(iRplot))])
        xlab{isub-1} = [' R = ' num2str(R(iRplot))];
    end
    %
    if 1 % wimax vs R,S, log, growthrate as colorscale
        hca = h(isub); isub=isub+1;             
        %[cs, hh] = contour(hca,R,S,x_imag_plot/normal_frequency,[-0.02:0.01:0.1]); clabel(cs, hh, 'labelspacing', 80);
        wiRSplot = squeeze(wimax(:,:,iTe2plot,iTiplot))*sqrt(1836);
        pcolor(hca,R,S,wiRSplot)
        %mesh(hca,R,S,x_imag_plot/normal_frequency)
        
        xlabel(hca,'R')
        ylabel(hca,'S')
        zlabel(hca,'\omega_{i}@\omega_{i,max}')
        ch = colorbar('peer',hca);
        %set(ch,'cscale','log')
        ylabel(ch,'\omega_{i}/\omega_{pi}@\omega_{i,max}')
        ylabel(ch,'\omega_{i}/\omega_{pi}')
        
        shading(hca,'flat')
        view(hca,[0 0 1])
        box(hca,'on')
        set(hca,'clim',5*[-1 1],'ylim',S([1 end]),'xlim',R([1 end]))
        set(ch,'ylim',[-2 5])
        title(hca,'Maximum growth rate')
    end
    if 1 % wimax vs R,S, log, growthrate as colorscale
        hca = h(isub); isub=isub+1;             
        %[cs, hh] = contour(hca,R,S,x_imag_plot/normal_frequency,[-0.02:0.01:0.1]); clabel(cs, hh, 'labelspacing', 80);
        wrRSplot = squeeze(wrmax(:,:,iTe2plot,iTiplot))*sqrt(1836);
        pcolor(hca,R,S,wrRSplot)
        %mesh(hca,R,S,x_imag_plot/normal_frequency)
        
        xlabel(hca,'R')
        ylabel(hca,'S')
        zlabel(hca,'\omega_{i}@\omega_{i,max}')
        ch = colorbar('peer',hca);
        %set(ch,'cscale','log')
        ylabel(ch,'\omega_{i}/\omega_{pi}@\omega_{i,max}')
        ylabel(ch,'\omega_{i}/\omega_{pi}')
        
        shading(hca,'flat')
        view(hca,[0 0 1])
        box(hca,'on')
        set(hca,'clim',10*[-1 1],'ylim',S([1 end]),'xlim',R([1 end]))
        set(ch,'ylim',[-2 5])
        title(hca,'Maximum growth rate')
    end
    if 1 % vph vs R,S, log, growthrate as colorscale
        hca = h(isub); isub=isub+1;     
        wavelengthRSplot = 2*pi./squeeze(kmax(:,:,iTe2plot,iTiplot));
        pcolor(hca,R,S,wavelengthRSplot)
        
        
            pcolor(hca,R,S,wavelengthRSplot)
            clim = [0 20];
            set(hca,'clim',clim,'ylim',S([1 end]),'xlim',R([1 end]),'zscale','lin')
            %set(hca,'zlim',[-0.01 0.12])
            ch = colorbar('peer',hca);                        
            ylabel(ch,'\lambda^*/\lambda_{De}')                
            set(ch,'ylim',clim)
                   
        shading(hca,'flat')
        view(hca,[0 0 1])
        box(hca,'on')
        grid(hca,'off')
        title(hca,'Wavelength at maximum growth rate')
        xlabel(hca,'R')
        ylabel(hca,'S')
        zlabel(hca,'\lambda^*')
        
    end
    if 1 % vph vs R,S, log, growthrate as colorscale
        hca = h(isub); isub=isub+1;     
        vphRSplot = squeeze(vphmax(:,:,iTe2plot,iTiplot));
        pcolor(hca,R,S,vphRSplot)
        
        
            pcolor(hca,R,S,vphRSplot)
            clim = [.1 1.2];
            set(hca,'clim',clim,'ylim',S([1 end]),'xlim',R([1 end]),'zscale','lin')
            %set(hca,'zlim',[-0.01 0.12])
            ch = colorbar('peer',hca);                        
            ylabel(ch,'v_{ph}/v_{te,bg}')                
            set(ch,'ylim',clim)
                   
        shading(hca,'flat')
        view(hca,[0 0 1])
        box(hca,'on')
        grid(hca,'off')
        title(hca,'Phase velocity at maximum growth rate')
        xlabel(hca,'R')
        ylabel(hca,'S')
        zlabel(hca,'\omega_{i,max}')
        
    end
    if 0 % vph vs R,S, log, growthrate as colorscale
        hca = h(isub); isub=isub+1;     
        surf2plot = x_imag_plot/normal_frequency;
        col2plot = vph_plot;
        col2plot(surf2plot<0) = NaN;
        doLog = 1;
        if doLog
            surf(hca,R,S,surf2plot,log10(col2plot))            
            clim = [-2 0];
            set(hca,'clim',clim,'ylim',S([1 end]),'xlim',R([1 end]),'zscale','lin')
            %set(hca,'zlim',[-0.01 0.12])
            ch = colorbar('peer',hca);                        
            ylabel(ch,'log_{10}(v_{ph}/v_{te,bg})')     
            %cticks = -flipdim([0 0.01:0.01:0.1 0.2:0.1:1 2],2);
            %set(ch,'ytick',cticks,'yscale','log')
            
        else
            surf(hca,R,S,surf2plot,col2plot)
            clim = [0 0.7];
            set(hca,'clim',clim,'ylim',S([1 end]),'xlim',R([1 end]),'zscale','lin')
            %set(hca,'zlim',[-0.01 0.12])
            ch = colorbar('peer',hca);                        
            ylabel(ch,'v_{ph}/v_{te,bg}')                
            set(ch,'ylim',clim)
        end            
        shading(hca,'flat')
        view(hca,[0 0 1])
        box(hca,'on')
        grid(hca,'off')
        title(hca,'Phase velocity at maximum growth rate')
        xlabel(hca,'R')
        ylabel(hca,'S')
        zlabel(hca,'\omega_{i,max}')
        
    end
    if 0 % e*Phi/kB*Tebg vs R,S
        hca = h(isub); isub=isub+1;     
        surf2plot = x_imag_plot/normal_frequency;
        col2plot = ephi;
        col2plot(surf2plot<0) = NaN;
        surf(hca,R,S,surf2plot,col2plot)
        %ctop=0.7;
        set(hca,'ylim',S([1 end]),'xlim',R([1 end]),'zscale','lin','clim',[0 0.8])
        %set(hca,'zlim',[-0.01 0.12])
        xlabel(hca,'R')
        ylabel(hca,'S')
        zlabel(hca,'\omega_{i,max}')
        ch = colorbar('peer',hca);
        %set(ch,'cscale','log')
        ylabel(ch,'v_{ph}/v_{te,bg}@\omega_{i,max}')
        ylabel(ch,'e\phi/k_BT_{e,bg}')
        shading(hca,'flat')
        view(hca,[0 0 1])
        box(hca,'on')
        grid(hca,'off')
        %set(ch,'ylim',[0 ctop])
        titlestr = {'Electrostatic potential required to capture beam e   ^-',...
            'e\phi = (v_{beam}-v_{te,beam}-v_{ph})^2 m_e/2'};
        title(hca,titlestr)
    end
        % add labels
    abc = {'a)','b)','c)','d)','e)','f)'};
    for kk = 1:nPlots        
        xl = get(h(kk),'xlim');
        yl = get(h(kk),'ylim');
        axes(h(kk))
        try abc{kk} = [abc{kk} xlab{kk}]; end
        text(xl(1)+diff(xl)*0.03,yl(end)-diff(xl)*0.03,abc{kk},'horizontalalignment','left','verticalalignment','top')
    end
    colormap(cn.cmap('bluered3'))
end

    