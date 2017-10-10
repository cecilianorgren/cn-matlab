% run_beam_solver

doLoad = 0;
doSave = 0;
savePath = '/Users/Cecilia/Research/EH2/Runs-solver_Buneman_kvR';
if doLoad
    dir(savePath)
    %load([savePath '/' 'Ti-2000eV_eventhickergrid'])
    %load('Ti-500eV_eventhickergrid')
else
    inputs = 2;
    switch inputs
        case 1 % more values 
            Te1 = 1600;
            Te2 = [120 60 20 10];
            Ti = [1600 3200 4800];
            S = 0.1:0.05:1.4;
            R = 0.09:0.05:0.99;
            k = 0.1:0.05:2; 
        case 2 % less values, ~10s 
            Te1 = 1600;
            Te2 = [20];
            Ti = [1600];
            S = [0.3 0.7 1.4];
            R = 0.09:0.2:0.99;
            k = 0.1:0.1:1.5;     
    end

    if 0 % rests from when considering external normalization 
        units = irf_units;
        e = units.e;
        eps0 = units.eps0;
        mp = units.mp;
        me = units.me;
        c = units.c;

        n = 0.06;
        omega_pe = sqrt(n*1e6*e^2/eps0/me); % rad/s
        vthe = c*sqrt(1-1./(Te1.*e./(me*c^2)+1).^2);
    end

    n = 0.06;
    [wr,wi] = beam_solver(Te1,Te2/Te1,Ti/Te1,R,S,k,n);
    [kmax,wrmax,wimax,vphmax,vph] = beam_solver_max(k,wr,wi,0);
    
    if doSave        
        saveName = [num2str(numel(k)) 'kx' num2str(numel(S)) 'Sx' num2str(numel(R)) 'Rx' num2str(numel(Te2)) 'Te2x' num2str(numel(Ti)) 'Ti'];
        save([savePath '/' saveName],'k','wr','wi','kmax','wrmax','wimax','vphmax','vph','Te1','Te2','Ti','R','S')
    end
end

%% calculate the additional quantities based on the beam model
[vphi,phi] = beam_solver_addons(vphmax,S,Te1);

%% plot results
% results stored in x_real_store(k,S,R) and x_imag_store(k,S,R)
% k vs R, loop over S
%% Surfaces
nPlots = 5;
nrows = 2;
for kk = 1:nPlots; h(kk)=subplot(nrows,ceil(nPlots/2),kk); end
for iR = 1
    isub = 1;
    if 1
        hca = h(isub); isub=isub+1;        
        surfc(hca,k,S,squeeze(x_real_store(:,:,iR))'/normal_frequency,squeeze(x_imag_store(:,:,iR))'/normal_frequency)
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'v_b/v_{te,bg}')
        zlabel(hca,'\omega_r/\omega_{pe}')
        ch = colorbar('peer',hca);
        set(hca,'clim',0.2*[-1 1])
        set(ch,'clim',[-1 1])
        title(hca,['R = ' num2str(R(iR)) ', T_{e,bg} = ' num2str(Te1) ' eV, T_{e,beam} = ' num2str(Te2) ' eV, T_i = ' num2str(Ti) ' eV'])
        ylabel(ch,'\omega_i/\omega_{pe}')
    end
    if 1
        hca = h(isub); isub=isub+1;        
        surfc(hca,k,S,squeeze(x_imag_store(:,:,iR))'/normal_frequency,squeeze(x_real_store(:,:,iR))'/normal_frequency)
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'v_b/v_{te,bg}')
        zlabel(hca,'\omega_i/\omega_{pe}')
        ch = colorbar('peer',hca);
        colormap(cn.cmap('islands'));
        ylabel(ch,'\omega_r/\omega_{pe}')
    end
    if 1
        hca = h(isub); isub=isub+1;        
        surfc(hca,k,S,squeeze(x_imag_store(:,:,iR))'/normal_frequency,squeeze(vph_store(:,:,iR))')
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'v_b/v_{te,bg}')
        zlabel(hca,'\omega_i/\omega_{pe}')
        ch = colorbar('peer',hca);
        colormap(cn.cmap('islands'));
        ylabel(ch,'v_{ph}/v_{te,bg}')
    end
    
    if 1
        hca = h(isub); isub=isub+1;        
        surfc(hca,k,S,squeeze(vph_store(:,:,iR))',squeeze(x_imag_store(:,:,iR))'/normal_frequency)
        xlabel(hca,'k\lambda_{De}')
        zlabel(hca,'v_{ph}/v_{te,bg}')
        ylabel(hca,'v_b/v_{te,bg}')
        ch = colorbar('peer',hca);
        colormap(cn.cmap('islands'));
        ylabel(ch,'\omega_i/\omega_{pe}')
        set(hca,'clim',0.2*[-1 1])
    end
    if 0
        hca = h(isub); isub=isub+1;        
        plot(hca,k,squeeze(x_real_store(:,:,iR))'/normal_frequency)
        xlabel(hca,'k\lambda_{De}')
        %ylabel(hca,'v_b/v_{te,bg}')
        %ch = colorbar('peer',hca);
        %set(hca,'clim',0.2*[-1 1])
        %set(ch,'clim',[-1 1])
    end
    if 0
        hca = h(isub); isub=isub+1;        
        plot(hca,k,squeeze(x_imag_store(:,:,iR))'/normal_frequency)
        xlabel(hca,'k\lambda_{De}')
        %ylabel(hca,'v_b/v_{te,bg}')
        %ch = colorbar('peer',hca);
    end
    
end

%% Lines, one plot for each R
units=irf_units;
doPrint=0;
nPlots = 4;
nrows = 2;
for kk = 1:nPlots; h(kk)=subplot(nrows,ceil(nPlots/2),kk); end
%set(gcf,'defaultcolororder',colors)
set(gcf,'position',[470 1 1201 955])
colors = jet(10);
set(0,'defaultAxesColorOrder',colors,...
      'defaultAxesLineStyleOrder','*|+|s')
colors = [0.4 0.4 1; 0.9 0 1;1 0 0;  1 0.7 0.3; .95 .95 0;0 0.8 0; 0 1 1];
%colors = spring(4);
set(0,'defaultAxesColorOrder',colors,...
      'defaultAxesLineStyleOrder','-|--|-.|:')
for iR = nR
    isub = 1;
    if 1 % wr
        hca = h(isub); isub=isub+1;        
        plot(hca,k,squeeze(x_real_store(:,:,iR))'/normal_frequency); hold(hca,'on');
        plot(hca,[0 k(end)],1*sqrt(units.me/units.mp)*[1 1],'k-');
        text(-0.23,sqrt(units.me/units.mp),'\omega_{pi}\rightarrow','parent',hca,'verticalalignment','middle')
        units=irf_units;        
        
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'\omega_r/\omega_{pe,tot}')  
        if max(get(hca,'ylim'))>2; ymaxx=2;
        else ymaxx = max(get(hca,'ylim'));
        end
        set(hca,'ylim',[0 ymaxx],'xlim',[0 2])
        title(hca,['R = ' num2str(R(iR),'%.2f')])
        legend(hca,[repmat('S=',numel(S),1) num2str(S','%.2f')],'location','bestoutside')
        hold(hca,'off')
    end
    if 1 % wi
        hca = h(isub); isub=isub+1; 
        plot(hca,k,squeeze(x_imag_store(:,:,iR))'/normal_frequency,'.'); hold(hca,'on')
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'\omega_i/\omega_{pe,tot}');
        plot(hca,[0 k(end)],0*[1 1],'k-');
        
        % indicate k*vd=omega_pe        
        k_B = 1./S/sqrt(2);
        plot(hca,k_B,0,'*');
        %text(2,c_s/vthe1,'\leftarrow c_{s,ie}','parent',hca,'verticalalignment','middle')

        
        wi_max = max(max(x_imag_store(:,:,iR)))/normal_frequency;
        if wi_max<0.1; wilim = ceil(wi_max*100)/100;
        else wilim = ceil(wi_max*10)/10; end
        if wilim > 1; wilim = 0.1; end
        %set(hca,'ylim',[0 max(get(hca,'ylim'))],'xlim',[0 2])
        set(hca,'ylim',[-0.01 wilim],'xlim',[0 2])
        legend(hca,[repmat('S=',numel(S),1) num2str(S','%.2f')],'location','bestoutside')
        hold(hca,'off')
    end
    if 1 % phase velocity
        hca=h(isub); isub=isub+1;
        v_prop=1;sqrt(2);
        plot(hca,k,v_prop*squeeze(vph_store(:,:,iR))','.'); hold(hca,'on')            
        ylabel(hca,'v_{ph}/v_{te,bg}','fontsize',14);
        xlabel(hca,'k\lambda_{De}','fontsize',14);  
        vlim = max(get(hca,'ylim')); vlim = 1.1*max(S); vlim = 1;
        if max(get(hca,'ylim'))>1; ymaxx=1;
        else ymaxx = max(get(hca,'ylim'));
        end        
        set(hca,'ylim',[0 vlim],'xlim',[0 2])
        
        legend(hca,[repmat('S=',numel(S),1) num2str(S','%.2f')],'location','bestoutside')
        
        plot(hca,[0 k(end)],vthi/vthe1*[1 1],'k-');
        text(0,vthi/vthe1,'v_{thi}\rightarrow','parent',hca,'verticalalignment','middle','horizontalalignment','right')
        
        units=irf_units;        
        % Buneman phase velocity at max growth for S = 1
        if 0
            vBS = 0.92;
            vBS = 2.65;
            v_B = (units.me/16/units.mp)^(1/3)*vBS;
            plot(hca,[0 k(end)],v_B*[1 1],'k-');
            text(0,v_B,['v_{B,S=' num2str(vBS) '} \rightarrow'],'parent',hca,'verticalalignment','middle','horizontalalignment','right')       
            vBS = 2.35;
        end
        k_B = 1./S/sqrt(2);
        for jj=1:nv;
            v_B = (units.me/16/units.mp)^(1/3)*S(jj);
            %plot(hca,[0 k(end)],v_B*[1 1],'k-');
            plot(hca,[k_B(jj) k_max_store(jj,iR)],[v_B vph_max_store(jj,iR)],'k-o');
            %text(min(k),v_B,['v_{B,S=' num2str(vBS) '} \rightarrow'],'parent',hca,'verticalalignment','middle','horizontalalignment','right')
        end
        title(hca,'o = v_{fluid,Buneman}')
        
        c_s = sqrt((3*Ti+Te2)*units.kB*11604/units.mp);
        plot(hca,[0 k(end)],c_s/vthe1*[1 1],'k-');
        text(1,c_s/vthe1,'\leftarrow c_{s,ie}','parent',hca,'verticalalignment','middle')

        c_s = sqrt(Te2*units.kB*11604/units.mp);
        plot(hca,[0 k(end)],c_s/vthe1*[1 1],'k-');        
        text(1,c_s/vthe1,'\leftarrow c_{s,e}','parent',hca,'verticalalignment','middle')
        
        plot(hca,k_max_store(:,iR),v_prop*vph_max_store(:,iR),'*k')
        text(2.1,vlim,'* v_{ph}@ \omega_{i,max}','parent',hca,'verticalalignment','bottom')
        
        

        hold(hca,'off')
    end
    if 1 % distribution
        hca = h(isub); isub=isub+1;  
        plot(hca,ves*1e-3,squeeze(fe_store(:,:,iR)),vis*1e-3,fi_store/10,'k-') 
        %[h,f,v] = whamp.plot_f(n(incl),mass(incl),t(incl),vd(incl),d(incl),a(incl),b(incl),pitchangle,plotoption);
        %legend(hca,num2str(S','%.1f'),'location','bestoutside')
        legend(hca,[repmat('S=',numel(S),1) num2str(S','%.2f')],'location','bestoutside')
        title(hca,['R = ' num2str(R(iR),'%.2f') ', T_{e,bg} = ' num2str(Te1) ' eV, T_{e,beam} = ' num2str(Te2) ' eV, T_i = ' num2str(Ti) ' eV'])
        xlabel(hca,'v [10^3 km/s]')
    end
    for ll=1:3; set(h(ll),'xlim',[min(k) max(k)]); end
    set(h(1),'ylim',[0 0.08])
    set(h(3),'ylim',[0.01 0.09])
    pause(1)
    if doPrint
        cn.print(['lines_Ti' num2str(Ti,'%.0f') 'eV__Te1' num2str(Te1,'%.0f') 'eV__Te2' num2str(Te2,'%.0f') 'eV__R' num2str(R(iR)*100,'%03.0f')])
    end
end

%% k vs S, loop over R
nPlots = 5;
nrows = 2;
for kk = 1:nPlots; h(kk)=subplot(nrows,ceil(nPlots/2),kk); end
for iv = 4
    wrmax = max(max(x_imag_store(:,iv,:)))/normal_frequency;
    isub = 1;
    if 1
        hca = h(isub); isub=isub+1;        
        surfc(hca,k,R,squeeze(x_real_store(:,iv,:))'/normal_frequency,squeeze(x_imag_store(:,iv,:))'/normal_frequency)
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'n_b/_{tot}')
        zlabel(hca,'\omega_r/\omega_{pe}')
        ch = colorbar('peer',hca);
        set(hca,'clim',1.2*wrmax*[-1 1])
        set(ch,'clim',[-1 1])
        title(hca,['S = ' num2str(S(iv)) ', T_{e,bg} = ' num2str(Te1) ' eV, T_{e,beam} = ' num2str(Te2) ' eV, T_i = ' num2str(Ti) ' eV'])
        ylabel(ch,'\omega_i/\omega_{pe}')
    end
    if 1
        hca = h(isub); isub=isub+1;        
        surfc(hca,k,R,squeeze(x_imag_store(:,iv,:))'/normal_frequency,squeeze(x_real_store(:,iv,:))'/normal_frequency)
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'n_b/_{tot}')
        zlabel(hca,'\omega_i/\omega_{pe}')
        ch = colorbar('peer',hca);
        colormap(cn.cmap('islands'));
        ylabel(ch,'\omega_r/\omega_{pe}')
    end
    if 1
        hca = h(isub); isub=isub+1;        
        surfc(hca,k,R,squeeze(x_imag_store(:,iv,:))'/normal_frequency,squeeze(vph_store(:,iv,:))')
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'n_b/_{tot}')
        zlabel(hca,'\omega_i/\omega_{pe}')
        ch = colorbar('peer',hca);
        colormap(cn.cmap('islands'));
        ylabel(ch,'v_{ph}/v_{te,bg}')
    end
    
    if 1
        hca = h(isub); isub=isub+1;        
        surfc(hca,k,R,squeeze(vph_store(:,iv,:))',squeeze(x_imag_store(:,iv,:))'/normal_frequency)
        xlabel(hca,'k\lambda_{De}')
        zlabel(hca,'v_{ph}/v_{te,bg}')
        ylabel(hca,'n_b/n_{tot}')
        ch = colorbar('peer',hca);
        colormap(cn.cmap('islands'));
        ylabel(ch,'\omega_i/\omega_{pe}')
        set(hca,'clim',1.2*wrmax*[-1 1])
    end
    if 0
        hca = h(isub); isub=isub+1;        
        plot(hca,k,squeeze(x_real_store(:,iv,:))'/normal_frequency)
        xlabel(hca,'k\lambda_{De}')
        %ylabel(hca,'v_b/v_{te,bg}')
        %ch = colorbar('peer',hca);
        %set(hca,'clim',0.2*[-1 1])
        %set(ch,'clim',[-1 1])
    end
    if 0
        hca = h(isub); isub=isub+1;        
        plot(hca,k,squeeze(x_imag_store(:,iv,:))'/normal_frequency)
        xlabel(hca,'k\lambda_{De}')
        %ylabel(hca,'v_b/v_{te,bg}')
        %ch = colorbar('peer',hca);
    end
    
end

%% Lines, one plot for each S
nPlots = 4;
nrows = 2;
for kk = 1:nPlots; h(kk)=subplot(nrows,ceil(nPlots/2),kk); end
set(gcf,'position',[470 1 1201 955])
colors = jet(4); %colors = colors(2:2:end,:);
colors = [0.4 0.4 1; 0.9 0 1;1 0 0;  1 0.7 0.3; .95 .95 0;0 0.8 0; 0 1 1];
%colors = spring(4);
set(0,'defaultAxesColorOrder',colors,...
      'defaultAxesLineStyleOrder','-|--|-.|:')
for iv = 1:nv
    isub = 1;
    if 1 % wr
        hca = h(isub); isub=isub+1;        
        plot(hca,k,squeeze(x_real_store(:,iv,:))'/normal_frequency)
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'\omega_r/\omega_{pe,tot}')
        wr_max = max (max(x_real_store(:,iv,:)))/normal_frequency;
        if wr_max>2; wrmaxlim=2;
        else wrmaxlim = ceil(wr_max*100)/100; end        
        set(hca,'ylim',[0 wrmaxlim],'xlim',[0 2])
        title(hca,['S = ' num2str(S(iv),'%.2f')])
        legend(hca,[repmat('R=',numel(R),1) num2str(R','%.2f')],'location','bestoutside')
    end
    if 1 % wi
        hca = h(isub); isub=isub+1;        
        plot(hca,k,squeeze(x_imag_store(:,iv,:))'/normal_frequency); hold(hca,'on')
        xlabel(hca,'k\lambda_{De}')
        ylabel(hca,'\omega_i/\omega_{pe,tot}')  
        plot(hca,[0 k(end)],0*[1 1],'k-'); hold(hca,'off')
        wi_max = max(max(x_imag_store(:,iv,:)))/normal_frequency;
        if wi_max < 0; wiminlim = floor(wi_max*100)/100;
        else wiminlim = -0.01; end
        if wi_max<0.1; wimaxlim = ceil(wi_max*100)/100;
        else wimaxlim = ceil(wi_max*10)/10; end
        if wi_max>1; wimaxlim = 0.1; end
        %set(hca,'ylim',[0 max(get(hca,'ylim'))],'xlim',[0 2])
        set(hca,'ylim',[wiminlim wimaxlim],'xlim',[0 2])
        legend(hca,[repmat('R=',numel(R),1) num2str(R','%.2f')],'location','bestoutside')
    end
    if 1 % phase velocity
        hca=h(isub); isub=isub+1;
        plot(hca,k,squeeze(vph_store(:,iv,:))'); hold(hca,'on')      
        ylabel(hca,'v_{ph}/v_{te,bg}','fontsize',14);
        xlabel(hca,'k\lambda_{De}','fontsize',14);          
        plot(hca,[0 k(end)],S(iv)*[1 1],'k-'); 
        text(0.25,S(iv),'v_{de}','parent',hca,'verticalalignment','bottom')
        plot(hca,[0 k(end)],(S(iv)*vthe1-vthe2)/vthe1*[1 1],'k-');
        text(0.25,(S(iv)*vthe1-vthe2)/vthe1,'v_{de}-v_{te,beam}','parent',hca,'verticalalignment','bottom')        
        plot(hca,[0 k(end)],vthi/vthe1*[1 1],'k-');
        text(-0.23,vthi/vthe1,'v_{thi}\rightarrow','parent',hca,'verticalalignment','middle')
        vlim = max(get(hca,'ylim')); vlim = 1.5; vlim = 1.1*S(iv);
        set(hca,'ylim',[0 vlim],'xlim',[0 2])
        plot(hca,k_max_store(iv,:),vph_max_store(iv,:),'*k')
        text(2.1,vlim,'* v_{ph}@ \omega_{i,max}','parent',hca,'verticalalignment','bottom')
        legend(hca,[repmat('R=',numel(R),1) num2str(R','%.2f')],'location','bestoutside')
        hold(hca,'off')
    end
    if 1 % distribution
        hca = h(isub); isub=isub+1;  
        plot(hca,ves,squeeze(fe_store(:,iv,:)),vis,fi_store/10,'k-'); hold(hca,'on')
        %plot(vph_max_store(iv,:))
        %v = logspace(3,5,100);
        %loglog(hca,v,cn.maxwellian(v,Te1,n*(1-R(1:2:end)))+cn.maxwellian(v,Te2,n*R))
        %[h,f,v] = whamp.plot_f(n(incl),mass(incl),t(incl),vd(incl),d(incl),a(incl),b(incl),pitchangle,plotoption);
        legend(hca,[repmat('R=',numel(R),1) num2str(R','%.2f')],'location','bestoutside')
        title(hca,['S = ' num2str(S(iv),'%.2f') ', T_{e,bg} = ' num2str(Te1) ' eV, T_{e,beam} = ' num2str(Te2) ' eV, T_i = ' num2str(Ti) ' eV'])
        
    end
    drawnow
    pause(2)
    %cn.print(['lines_Ti-' num2str(Ti,'%.0f') 'eV__Te1-' num2str(Te1,'%.0f') 'eV_Te2-' num2str(Te2,'%.0f') 'eV__S' num2str(S(iv)*100,'%03.0f')])
end

%% Plot a surface plot of vph@maxwi vs R and S
vph_plot=vph_max_store;
vph_plot(vph_max_store>2)=NaN;
vph_plot(vph_max_store<1e-3)=NaN;
x_imag_plot = x_imag_max_store_2Dvalues;
x_imag_plot(x_imag_max_store_2Dvalues/normal_frequency>0.2)=NaN;
x_imag_plot(vph_max_store>2)=NaN;

nPlots = 6;
nrows = 2;
for kk = 1:nPlots; h(kk)=subplot(nrows,ceil(nPlots/2),kk); end
isub=1;
if 1 % vph vs R,S, lin
    hca = h(isub); isub=isub+1;         
    surfc(hca,R,S,vph_plot)
    set(hca,'clim',[0 1],'ylim',S([1 end]),'xlim',R([1 end]),'zlim',[-0.01 1])
    xlabel(hca,'R')
    ylabel(hca,'S')
    zlabel(hca,'v_{ph}@\omega_{i,max}')
    ch = colorbar('peer',hca);
    ylabel(ch,'v_{ph}/v_{te,bg}')
end
if 1 % vph vs R,S, log
    hca = h(isub); isub=isub+1;     
    surfc(hca,R,S,vph_plot)
    set(hca,'clim',[0 1],'ylim',S([1 end]),'xlim',R([1 end]),'zscale','log')
    xlabel(hca,'R')
    ylabel(hca,'S')
    zlabel(hca,'v_{ph}@\omega_{i,max}')
    ch = colorbar('peer',hca);
    %set(ch,'cscale','log')
    ylabel(ch,'v_{ph}/v_{te,bg}')
end
if 1 % vph vs R,S, log
    hca = h(isub); isub=isub+1;     
    surfc(hca,R,S,log(vph_plot))
    set(hca,'ylim',S([1 end]),'xlim',R([1 end]),'zscale','lin')
    xlabel(hca,'R')
    ylabel(hca,'S')
    zlabel(hca,'log(v_{ph}@\omega_{i,max})')
    ch = colorbar('peer',hca);
    ylabel(ch,'v_{ph}/v_{te,bg}')    
end
if 1 % vph vs R,S, log, growthrate as colorscale
    hca = h(isub); isub=isub+1;     
    surfc(hca,R,S,vph_plot,x_imag_max_store_2Dvalues/normal_frequency)
    set(hca,'clim',[0 0.1],'ylim',S([1 end]),'xlim',R([1 end]),'zscale','log')
    xlabel(hca,'R')
    ylabel(hca,'S')
    zlabel(hca,'v_{ph}@\omega_{i,max}')
    ch = colorbar('peer',hca);
    %set(ch,'cscale','log')
    ylabel(ch,'\omega_{i}@\omega_{i,max}')
end

if 1 % vph vs R,S, log, growthrate as colorscale
    hca = h(isub); isub=isub+1;     
    surfc(hca,R,S,x_imag_plot/normal_frequency)
    set(hca,'clim',[0 0.1],'ylim',S([1 end]),'xlim',R([1 end]),'zlim',[-0.01 0.12],'zscale','lin')
    xlabel(hca,'R')
    ylabel(hca,'S')
    zlabel(hca,'\omega_{i}@\omega_{i,max}')
    ch = colorbar('peer',hca);
    %set(ch,'cscale','log')
    ylabel(ch,'\omega_{i}@\omega_{i,max}')
    shading(hca,'flat')
end

if 1 % vph vs R,S, log, growthrate as colorscale
    hca = h(isub); isub=isub+1;     
    surf2plot = x_imag_plot/normal_frequency;
    col2plot = vph_plot;
    col2plot(surf2plot<0) = NaN;
    surfc(hca,R,S,surf2plot,col2plot)
    set(hca,'clim',[0 1],'ylim',S([1 end]),'xlim',R([1 end]),'zlim',[-0.01 0.12],'zscale','lin')
    xlabel(hca,'R')
    ylabel(hca,'S')
    zlabel(hca,'\omega_{i,max}')
    ch = colorbar('peer',hca);
    %set(ch,'cscale','log')
    ylabel(ch,'v_{ph}@\omega_{i,max}')
end
% add tielt to first panel
title(h(1),['T_{e,bg} = ' num2str(Te1) ' eV, T_{e,beam} = ' num2str(Te2) ' eV, T_i = ' num2str(Ti) ' eV'])

%% Prepare xes1 input file
units = irf_units;
RR = [0.1 0.15 0.2 0.2 0.2 0.2 0.4 0.5  0.5  0.5 0.6 0.98]; 
SS = [0.4 0.7 0.7 0.4 1.5 0.5 0.5 0.55 0.65 0.8 0.6 0.9 ];
disp('----')
for ir = 1:numel(RR)    
    ope1 = omega_pe*sqrt(1-RR(ir));
    ope2 = omega_pe*sqrt(RR(ir));
    opi = omega_pe*sqrt(units.me/units.mp);
    vte1 = irf_plasma_calc(B,n,no,Te1,Ti,'Vte'); % m/s 
    vte2 = irf_plasma_calc(B,n,no,20,Ti,'Vte'); % m/s 
    vti = irf_plasma_calc(B,n,no,Te1,Ti,'Vtp'); % m/s 
    vde2 = SS(ir)*vte1;

    xes1_ope1 = ope1/omega_pe;
    xes1_ope2 = ope2/omega_pe;
    xes1_opi = opi/omega_pe;
    xes1_vte1 = vte1/vte1;
    xes1_vte2 = vte2/vte1;
    xes1_vti = vti/vte1;
    xes1_vde2 = vde2/vte1;
    
    disp(['R = ' num2str(RR(ir),'%.2f') ', S = ' num2str(SS(ir),'%.2f') ', wpe1 = ' num2str(xes1_ope1,'%.3f') ', wpe2 = ' num2str(xes1_ope2,'%.3f') ', wpi = ' num2str(xes1_opi,'%.3f') ', vte1 = ' num2str(xes1_vte1,'%.3f') ', vte2 = ' num2str(xes1_vte2,'%.3f') ', vti = ' num2str(xes1_vti,'%.3f') ', vde2 = ' num2str(xes1_vde2,'%.3f')])
end




