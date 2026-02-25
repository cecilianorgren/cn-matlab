cd /Users/Cecilia/Research/EH2/Runs-solver_Buneman_kvR/
doLoad = 1;
doClearup = 1;
doPlot = 1;
doPhi = 1;
if doLoad
    load('Ti-2000eV_eventhickergrid')
    %load('Ti-500eV_eventhickergrid')
    nR = numel(R);
    nv = numel(S);
    nk = numel(k)

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
if doClearup
    %% clear up the matrix from phoney values
    k_store = repmat(k',1,nv,nR);
    vph_store = x_real_store/normal_frequency./k_store/sqrt(2);
    x_imag_store(x_real_store<0) = NaN;
    vph_store(x_real_store<0) = NaN;
    x_real_store(x_real_store<0) = NaN;
    x_real_store(x_real_store/normal_frequency>2.4) = NaN;
    %% find the point where the growth is max
    % find the max for each R and S, i.e. one per k-vector
    x_imag_max_store_3D = zeros(size(x_imag_store));
    x_imag_max_store_2D = zeros(nv,nR);
    x_imag_max_store_2Dvalues = zeros(nv,nR);
    vph_max_store = zeros(nv,nR);
    k_max_store  = zeros(nv,nR);
    %(x_imag_max_store_2D(iv,:)
    for iR = 1:nR
        for iS = 1:nv
            [ind,~] = find(x_imag_store(:,iS,iR)==max(x_imag_store(:,iS,iR)));

            x_imag_max_store_3D(ind,iS,iR)=1;    
            if ~isempty(ind)
                x_imag_max_store_2D(iS,iR)=ind;
                vph_max_store(iS,iR) = vph_store(ind,iS,iR);
                k_max_store(iS,iR) = k_store(ind,iS,iR);
                x_imag_max_store_2Dvalues(iS,iR) = x_imag_store(ind,iS,iR);
            else
                x_imag_max_store_2D(iS,iR)=NaN;
                vph_max_store(iS,iR) = NaN;
                k_max_store(iS,iR) = NaN;
                x_imag_max_store_2Dvalues(iS,iR) = NaN;
            end
            disp(['iR=' num2str(iR) ', iS=' num2str(iS) ', max_ind=' num2str(ind) ', vph_max_ind=' num2str(vph_max_store(iS,iR))])        
        end
    end
end
if doPlot
    vph_plot=vph_max_store;
    vph_plot(vph_max_store>2)=NaN;
    vph_plot(vph_max_store<1e-3)=NaN;
    x_imag_plot = x_imag_max_store_2Dvalues;
    x_imag_plot(x_imag_max_store_2Dvalues/normal_frequency>0.2)=NaN;
    x_imag_plot(vph_max_store>2)=NaN;
    x_imag_plot(x_imag_max_store_2Dvalues/normal_frequency<-0.02)=NaN;
    
    if doPhi
        %%
        units = irf_units;
        vthebeam = cn_eV2v(Te2,'eV'); % km/s
        vthebg = cn_eV2v(Te1,'eV'); % km/s
        vbeam = repmat(S',1,numel(R))*vthebg;
        ephi = (vbeam-1*vthebeam-vph_plot*vthebg).^2*1e6*units.me/2/units.e/Te1;
    end
%%    
    set(0,'defaultAxesFontSize',16);
    set(0,'DefaultTextFontSize',16);
    set(0,'defaultAxesFontUnits','pixels');
    set(0,'defaultTextFontUnits','pixels');
    nPlots = 3;
    nrows = 1;
    for kk = 1:nPlots; h(kk)=subplot(nrows,ceil(nPlots/2+1),kk); end
    isub=1;
    if 1 % vph vs R,S, log, growthrate as colorscale
        hca = h(isub); isub=isub+1;             
        %[cs, hh] = contour(hca,R,S,x_imag_plot/normal_frequency,[-0.02:0.01:0.1]); clabel(cs, hh, 'labelspacing', 80);
        pcolor(hca,R,S,x_imag_plot/normal_frequency)
        %mesh(hca,R,S,x_imag_plot/normal_frequency)
        
        xlabel(hca,'R')
        ylabel(hca,'S')
        zlabel(hca,'\omega_{i}@\omega_{i,max}')
        ch = colorbar('peer',hca);
        %set(ch,'cscale','log')
        ylabel(ch,'\omega_{i}/\omega_{pe}@\omega_{i,max}')
        ylabel(ch,'\omega_{i}/\omega_{pe}')
        
        shading(hca,'flat')
        view(hca,[0 0 1])
        box(hca,'on')
        set(hca,'clim',0.1*[-1 1],'ylim',S([1 end]),'xlim',R([1 end]),'zlim',[-0.01 0.12],'zscale','lin')
        set(ch,'ylim',[-0.03 0.1])
        title(hca,'Maximum growth rate')
    end
    if 1 % vph vs R,S, log, growthrate as colorscale
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
    if 1 % e*Phi/kB*Tebg vs R,S
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
    colormap(cn.cmap('bluered3'))
end

    