cd /Users/Cecilia/Research/EH2/Runs-solver_Buneman_kvR/
doLoad = 1;
doClearup = 1;

if doLoad
    %load('Ti-2000eV_eventhickergrid')
    %load('Ti-500eV_eventhickergrid')
    load('Te1-2000_Te2-120_Ti-2000eV_R15S15k40')
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
vph_plot = vph_max_store;
vph_plot(vph_max_store>2) = NaN;
vph_plot(vph_max_store<1e-3) = NaN;


f_tot = 1;
f_vthebeam = 1;
units = irf_units;
vthebeam = cn_eV2v(Te2,'eV'); % km/s
vthebg = cn_eV2v(Te1,'eV'); % km/s
vbeam = repmat(S',1,numel(R))*vthebg;

vphi = f_tot*(vbeam-f_vthebeam*vthebeam-vph_plot*vthebg);
% Trapping velocity from Omura1996
Rmat = repmat(R,numel(S),1);

a1 = vbeam-1*vthebeam;
a2 = 8/3*Rmat;
a3 = 9*vbeam*vthebeam/2./Rmat./a1./a1;
vT = a1*sqrt(a2*(sqrt(1+a3)-1))/1000;
%vT = (vbeam-1*vthebeam)*(sqrt(8/3*Rmat*sqrt(1+9*vbeam*vthebeam/2./Rmat/(vbeam-vthebeam)/(vbeam-vthebeam))-1))
% vphi = vT;
vph = vph_plot*vthebg;

titlestr = ['T_{e,bg} = ' num2str(Te1) ' eV, T_{e,beam} = ' num2str(Te2) ' eV, T_i = ' num2str(Ti) ' eV'];
titlestr2 = ['v_{\phi} = ' num2str(f_tot,'%.1f') '\times (v_{beam} - ' num2str(f_vthebeam,'%.1f') '*v_{the,beam} - v_{ph})'];
%%

% daniels fit
pexp=1;
pcoeff=15;
dppxx = linspace(1e-4,1e2,10);
dppyy = pcoeff*dppxx.^pexp;
dpstr = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];

% fit to model
pexp = 1;
pcoeff = 3;
    

% set up plot
nPlots = 3;
nrows = 1;
ncols = ceil(nPlots/nrows);
for kk = 1:nPlots; h(kk)=subplot(1,ncols,kk); end
isub = 1;

if 0
    hca = h(isub); isub = isub+1;
    ppxx = linspace(1e-3,1e1,1000);
    pexp = -1;
    pcoeff = 5e-1;
    ppyy = pcoeff*(ppxx).^(pexp);
    pstr = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];
    plx1 = size(vph,2); ply1 = 1:size(vph,1); % max density
    plx2 = 1; ply2 = 1:size(vph,1); % min density
    loglog(hca,ppxx,ppyy,'k',...
           vph(plx1,ply1)./vthebg,vphi(plx1,ply1)./vph(plx1,ply1),'r--',...
           vph(plx2,ply2)./vthebg,vphi(plx2,ply2)./vph(plx2,ply2),'r--',...
           vph./vthebg,vphi./vph,'o')
    set(hca,'xlim',[1e-3 1e1],'ylim',[1e-2 1e2])
    xlabel(hca,'v_{ph}/v_{te,bg}')
    ylabel(hca,'v_{\phi}/v_{ph}')
    legend(hca,pstr,['R=' num2str(R(end))],'location','best')
    title(hca,titlestr)
end
if 0
    hca = h(isub); isub = isub+1;
    
    xlim = [1e-3 1e1];
    ylim = [1e-3 1e2];
    ppxx = linspace(1e-4,1e1,1000);
    pexp = 1;
    pcoeff = 0.5e1;
    ppyy = pcoeff*(ppxx).^(pexp);
    pstr = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];

    plx1 = size(vph,2); ply1 = 1:size(vph,1); % max density
    plx2 = 1; ply2 = 1:size(vph,1); % min density          

    loglog(ppxx,ppyy,'k',...
           dppxx,dppyy,'g-',...
           vph(plx1,ply1)./vthebg,vph(plx1,ply1)./vphi(plx1,ply1),'r--',...
           vph(plx2,ply2)./vthebg,vph(plx2,ply2)./vphi(plx2,ply2),'r--',...
           vph(plxs,plys)./vthebg,vph(plxs,plys)./vphi(plxs,plys),'o')
    legend(pstr,dpstr,['R=' num2str(R(end))],'location','best')

    set(gca,'xlim',xlim,'ylim',ylim)
    xlabel('v_{ph}/v_{te,bg}')
    ylabel('v_{ph}/v_{\phi}')

    title(titlestr)
end
if 1
    hca = h(isub); isub = isub+1;
    xlim = [1e-3 1e1];
    ylim = [1e-3 1e2];
    ppxx = linspace(1e-4,1e1,1000);
    %pexp = 1;
    %pcoeff = 1.2e1;
    ppyy = pcoeff*(ppxx).^(pexp);
    pstr = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];


    plx1 = size(vph,2); ply1 = 1:size(vph,1); % max density
    plx2 = 1; ply2 = 1:size(vph,1); % min density
    plxs = 1:size(vph,2); plys = 1:size(vph,1); % all R and S
    plys = 1:floor(size(vph,2)/3); plxs = 1:size(vph,1); % small Rs
    plym = ceil(size(vph,2)/3):floor(size(vph,2)*2/3); plxm = 1:size(vph,1); % middle Rs
    plyh = ceil(size(vph,2)*2/3):size(vph,2); plxh = 1:size(vph,1); % big Rs
   
    loglog(hca,ppxx,ppyy,'k',...
           dppxx,dppyy,'g-',...
           ...%vph(plx1,ply1)./vthebg,vph(plx1,ply1)./vphi(plx1,ply1),'b--',...
           ...%vph(plx2,ply2)./vthebg,vph(plx2,ply2)./vphi(plx2,ply2),'k--',...
           vph(plxs,plys)./vthebg,vph(plxs,plys)./vphi(plxs,plys),'ko',...
           vph(plxm,plym)./vthebg,vph(plxm,plym)./vphi(plxm,plym),'ro',...
           vph(plxh,plyh)./vthebg,vph(plxh,plyh)./vphi(plxh,plyh),'bo')
    Rstrs = ['R=' num2str(R(plys(1)),'%.2f') '-' num2str(R(plys(end)),'%.2f')];
    Rstrm = ['R=' num2str(R(plym(1)),'%.2f') '-' num2str(R(plym(end)),'%.2f')];
    Rstrh = ['R=' num2str(R(plyh(1)),'%.2f') '-' num2str(R(plyh(end)),'%.2f')];
       
    set(hca,'colororder',[0 0 0;1 0 0; 0 0 1])
    irf_legend(hca,{Rstrs,Rstrm,Rstrh},[0.95 0.05])
    legend(hca,pstr,dpstr,'location','northwest')
    
    %text(xlim(end),ylim(1),{Rstrs,Rstrm,Rstrh},'horizontalalignment','right','verticalalignment','bottom')
    
    set(hca,'xlim',xlim,'ylim',ylim)
    xlabel(hca,'v_{ph}/v_{te,bg}')
    ylabel(hca,'v_{ph}/v_{\phi}')

    title(hca,titlestr)
end
if 1
    hca = h(isub); isub = isub+1;
    xlim = [1e-3 1e1];
    ylim = [1e-3 1e2];
    ppxx = linspace(1e-4,1e1,1000);
    %pexp = 1;
    %pcoeff = 1.2e1;
    ppyy = pcoeff*(ppxx).^(pexp);
    pstr = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];


    plx1 = size(vph,2); ply1 = 1:size(vph,1); % max density
    plx2 = 1; ply2 = 1:size(vph,1); % min density
    plxs = 1:size(vph,2); plys = 1:size(vph,1); % all R and S
    plxs = 1:floor(size(vph,2)/3); plys = 1:size(vph,1); % small Rs
    plxm = ceil(size(vph,2)/3):floor(size(vph,2)*2/3); plym = 1:size(vph,1); % middle Rs
    plxh = ceil(size(vph,2)*2/3):size(vph,2); plyh = 1:size(vph,1); % big Rs


    loglog(hca,ppxx,ppyy,'k',...
           dppxx,dppyy,'g-',...
           ...%vph(plx1,ply1)./vthebg,vph(plx1,ply1)./vphi(plx1,ply1),'b--',...
           ...%vph(plx2,ply2)./vthebg,vph(plx2,ply2)./vphi(plx2,ply2),'k--',...
           vph(plxs,plys)./vthebg,vph(plxs,plys)./vphi(plxs,plys),'ko',...
           vph(plxm,plym)./vthebg,vph(plxm,plym)./vphi(plxm,plym),'ro',...
           vph(plxh,plyh)./vthebg,vph(plxh,plyh)./vphi(plxh,plyh),'bo')
     
    Sstrs = ['S=' num2str(S(plxs(1)),'%.2f') '-' num2str(S(plxs(end)),'%.2f')];
    Sstrm = ['S=' num2str(S(plxm(1)),'%.2f') '-' num2str(S(plxm(end)),'%.2f')];
    Sstrh = ['S=' num2str(S(plxh(1)),'%.2f') '-' num2str(S(plxh(end)),'%.2f')];
    
    set(hca,'colororder',[0 0 0;1 0 0; 0 0 1])
    irf_legend(hca,{Sstrs,Sstrm,Sstrh},[0.95 0.05])
    legend(hca,pstr,dpstr,'location','best')
    
    %text(xlim(end),ylim(1),{Rstrs,Rstrm,Rstrh},'horizontalalignment','right','verticalalignment','bottom')
    
    set(hca,'xlim',xlim,'ylim',ylim)
    xlabel(hca,'v_{ph}/v_{te,bg}')
    ylabel(hca,'v_{ph}/v_{\phi}')

    title(hca,titlestr2)
end
if 1
    hca = h(isub); isub = isub+1;
    xlim = [1e-3 1e1];
    ylim = [1e-3 1e2];
    ppxx = linspace(1e-4,1e1,1000);
    %pexp = 1;
    %pcoeff = 0.5e1;
    ppyy = pcoeff*(ppxx).^(pexp);
    pstr = ['fit: y = ' num2str(pcoeff) '*x^{' num2str(pexp) '}'];

    plxs = 1:size(vph,2); plys = 1:size(vph,1); % all R and S
    step=2;
    plys = 1:step:floor(size(vph,2)/3); 
    plxs = 1:step:floor(size(vph,2)/3); % small Rs
    plym = ceil(size(vph,2)/3):step:floor(size(vph,2)*2/3); 
    plxm = ceil(size(vph,2)/3):step:floor(size(vph,2)*2/3);  % middle Rs
    plyh = ceil(size(vph,2)*2/3):step:size(vph,2); 
    plxh = ceil(size(vph,2)*2/3):step:size(vph,2); % big Rs
   
    if 0
        loglog(hca,ppxx,ppyy,'k',...
               dppxx,dppyy,'g-',...
               vph(plxs,plys)./vthebg,vph(plxs,plys)./vphi(plxs,plys),'ko',...
               vph(plxm,plys)./vthebg,vph(plxm,plys)./vphi(plxm,plys),'ro',...
               vph(plxh,plys)./vthebg,vph(plxh,plys)./vphi(plxh,plys),'bo',...
               vph(plxs,plym)./vthebg,vph(plxs,plym)./vphi(plxs,plym),'c.',...
               vph(plxm,plym)./vthebg,vph(plxm,plym)./vphi(plxm,plym),'y.',...
               vph(plxh,plym)./vthebg,vph(plxh,plym)./vphi(plxh,plym),'g.',...
               vph(plxs,plyh)./vthebg,vph(plxs,plyh)./vphi(plxs,plyh),'kx',...
               vph(plxm,plyh)./vthebg,vph(plxm,plyh)./vphi(plxm,plyh),'rx',...
               vph(plxh,plyh)./vthebg,vph(plxh,plyh)./vphi(plxh,plyh),'bx')

        Sstrs = ['S=' num2str(S(plxs(1)),'%.2f') '-' num2str(S(plxs(end)),'%.2f')];
        Sstrm = ['S=' num2str(S(plxm(1)),'%.2f') '-' num2str(S(plxm(end)),'%.2f')];
        Sstrh = ['S=' num2str(S(plxh(1)),'%.2f') '-' num2str(S(plxh(end)),'%.2f')];    

        set(hca,'colororder',[0 0 0;1 0 0; 0 0 1])
        irf_legend(hca,{Sstrs,Sstrm,Sstrh},[0.95 0.05])
        legend(hca,pstr,dpstr,'location','best')
    else
        loglog(hca,ppxx,ppyy,'k',...
               dppxx,dppyy,'g-',...
               ...%vph(plx1,ply1)./vthebg,vph(plx1,ply1)./vphi(plx1,ply1),'b--',...
               ...%vph(plx2,ply2)./vthebg,vph(plx2,ply2)./vphi(plx2,ply2),'k--',...
               vph(plxs,plys)./vthebg,vph(plxs,plys)./vphi(plxs,plys),'ks',...
               vph(plxm,plys)./vthebg,vph(plxm,plys)./vphi(plxm,plys),'rs',...
               vph(plxh,plys)./vthebg,vph(plxh,plys)./vphi(plxh,plys),'bs',...
               vph(plxs,plym)./vthebg,vph(plxs,plym)./vphi(plxs,plym),'k.',...
               vph(plxm,plym)./vthebg,vph(plxm,plym)./vphi(plxm,plym),'r.',...
               vph(plxh,plym)./vthebg,vph(plxh,plym)./vphi(plxh,plym),'b.',...
               vph(plxs,plyh)./vthebg,vph(plxs,plyh)./vphi(plxs,plyh),'ko',...
               vph(plxm,plyh)./vthebg,vph(plxm,plyh)./vphi(plxm,plyh),'ro',...
               vph(plxh,plyh)./vthebg,vph(plxh,plyh)./vphi(plxh,plyh),'bo')

        Sstrs = ['S=' num2str(S(plxs(1)),'%.2f') '-' num2str(S(plxs(end)),'%.2f')];
        Sstrm = ['S=' num2str(S(plxm(1)),'%.2f') '-' num2str(S(plxm(end)),'%.2f')];
        Sstrh = ['S=' num2str(S(plxh(1)),'%.2f') '-' num2str(S(plxh(end)),'%.2f')];    
        Rstrs = ['R=' num2str(R(plxs(1)),'%.2f') '-' num2str(R(plxs(end)),'%.2f')];
        Rstrm = ['R=' num2str(R(plxm(1)),'%.2f') '-' num2str(R(plxm(end)),'%.2f')];
        Rstrh = ['R=' num2str(R(plxh(1)),'%.2f') '-' num2str(R(plxh(end)),'%.2f')];
    
        set(hca,'colororder',[0 0 0;1 0 0; 0 0 1])
        irf_legend(hca,{Sstrs,Sstrm,Sstrh},[0.95 0.05])
        legend(hca,pstr,dpstr,'location','best')
        %text(xlim(1)*1.1,ylim(end)*0.9,{[Rstrs],[Rstrm],[Rstrh]},'horizontalalignment','left','verticalalignment','top')
        if 0
            legind = 1;
            for ii = 1:3;
                for jj = 1:3;
                    legstr{legind} = ['R=' num2str(R(plys(jj)),'%.2f') '-' num2str(R(plys(end)),'%.2f') ', S=' num2str(S(plys(1)),'%.2f') '-' num2str(S(plys(end)),'%.2f')];
                    legind=legind+1;
                end
            end
        end
    end    
        
        
        

              
    
    
    %text(xlim(end),ylim(1),{Rstrs,Rstrm,Rstrh},'horizontalalignment','right','verticalalignment','bottom')
    
    set(hca,'xlim',xlim,'ylim',ylim)
    xlabel(hca,'v_{ph}/v_{te,bg}')
    ylabel(hca,'v_{ph}/v_{\phi}')

    title(hca,'Colors: S, symbols: R')
end

%% make fit
fitx = reshape(vph./vthebg,numel(vph./vthebg),1);
fity = reshape(vph./vphi,numel(vph./vphi),1);
nanind= isnan(fitx);
fitt = fit(fitx(~nanind),fity(~nanind),'power1');
hold(gca,'on')
ppyy = fitt.a*ppxx.^fitt.b;
loglog(ppxx,ppyy)
hold(gca,'off')



