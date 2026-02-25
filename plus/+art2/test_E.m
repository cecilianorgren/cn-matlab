% run art2.E

% Model
event = 15; % eh, sc
switch event
    case 13
        ic = 3;
        peroffs = 3;
        paroffs = 0; 
        phi0 = 200;        
        lr = 4;
        lz = 5;
        r = 2;
        t0 = toepoch([2007 08 31 10 17 51.794]);
        tint = t0 + [-0.03 0.03];
        v = 436;
        R=0;
    case 14    
        ic = 4;
        peroffs = 1;
        paroffs = -2;        
        phi0 = 210;
        lr = 9;
        lz = 7;
        r =3.6; % distance from center of the eh
        t0 = toepoch([2007 08 31 10 17 51.728]);
        tint = t0 + [-0.04 0.04];
        v = 436;
        R=0;
    case 15
        ic = 3;
        peroffs = 2;
        paroffs = 0; 
        phi0 = 200;        
        lr = 8;
        lz = 5;
        r = 6.3;
        t0 = toepoch([2007 08 31 10 17 51.794]);
        tint = t0 + [-0.03 0.03];
        v = 436;
        R=0.1;
    case 23
        ic = 3;
        peroffs = 0;
        paroffs = 4;        
        phi0 = 230;
        lr = 9;
        lz = 7;
        r = 1; % distance from center of the eh
        t0 = toepoch([2007 08 31 10 17 45.738]);
        tint = t0 + [-0.08 0.08];
        v = 436;
        R=0;
    case 24
        ic = 4;
        peroffs = 1;
        paroffs = 1;        
        phi0 = 210;
        lr = 9;
        lz =7;
        r =1; % distance from center of the eh
        t0 = toepoch([2007 08 31 10 17 45.662]);
        tint = t0 + [-0.04 0.04];
        v = 436;
        R=0;
    case 34
        ic = 4;
        peroffs = 1;
        paroffs = 0;        
        phi0 = 130;
        lr = 9;
        lz =4;
        r =6; % distance from center of the eh
        t0 = toepoch([2007 08 31 10 17 46.494]);
        tint = t0 + [-0.04 0.04];
        v = 436;  
        R = 1;
end
nr = 1;
zl = 3*lz; nz = 100; z = linspace(-zl,zl,nz);

for ii=1:nr
    [Ertmp Eztmp] = art2.E(phi0,r,-z,lr,lz);
    Er{ii} = Ertmp; Ez{ii} = Eztmp;
end
clear Ertmp Eztmp;


if 0
    % Figure
    figure(61)
    set(gcf,'position',[560 322 463 626]);
    h(1) = subplot(2,1,1); h(2) = subplot(2,1,2);
    hca = h(1); plot(hca,z,Ez{1}+Er{1}); xlabel(hca,'z [km]'); ylabel(hca,'E_z'); set(hca,'xlim',[-zl zl]); hold(hca,'on');
    hca = h(2); plot(hca,z,Er{1}); xlabel(hca,'z [km]'); ylabel(hca,'E_r'); set(hca,'xlim',[-zl zl]); hold(hca,'on');

    % Real field
    c_eval('Epar = irf_tlim(Epar?,tint);',ic);
    c_eval('Eper = irf_tlim(Eper?,tint);',ic);
    zE = Epar(:,1)*v; zE = zE - mean(zE);

    hca = h(1); plot(hca,zE,Epar(:,2)-paroffs,'r'); 
    hca = h(2); plot(hca,zE,Eper(:,2)-peroffs,'r'); 

    hold(h(1),'off'); hold(h(2),'off');

    strTitle={['l_{\perp} = ' num2str(lr) ' km,  l_{||} = ' num2str(lz) ' km,  crossing at r =' num2str(r) ' km'],...
              ['v = ' num2str(v) ' km/s,  t_0 = ' irf_time(t0,'isoshort')],...
              ['C' num2str(ic) ',  \phi_0 = ' num2str(phi0) ' V,  offs_{\perp} = ' num2str(peroffs) ' mV/m,  offs_{||} = ' num2str(paroffs) ' mV/m']};
    title(h(1),strTitle)
else
    % Add a part of Eperp to Epar, error coming in while constructing Epar from Esp    
    % Figure
    figure(61)
    set(gcf,'position',[560 322 463 626]);
    h(1) = subplot(2,1,1); h(2) = subplot(2,1,2);
    hca = h(1); plot(hca,z,Ez{1}+R*Er{1}); xlabel(hca,'z [km]'); ylabel(hca,'E_z'); set(hca,'xlim',[-zl zl]); hold(hca,'on');
    hca = h(2); plot(hca,z,Er{1}); xlabel(hca,'z [km]'); ylabel(hca,'E_r'); set(hca,'xlim',[-zl zl]); hold(hca,'on');

    % Real field
    c_eval('Epar = irf_tlim(Epar?,tint);',ic);
    c_eval('Eper = irf_tlim(Eper?,tint);',ic);
    zE = Epar(:,1)*v; zE = zE - mean(zE);

    hca = h(1); plot(hca,zE,Epar(:,2)-paroffs,'r'); 
    hca = h(2); plot(hca,zE,Eper(:,2)-peroffs,'r'); 

    hold(h(1),'off'); hold(h(2),'off');

    strTitle={['l_{\perp} = ' num2str(lr) ' km,  l_{||} = ' num2str(lz) ' km,  crossing at r =' num2str(r) ' km'],...
              ['v = ' num2str(v) ' km/s,  t_0 = ' irf_time(t0,'isoshort')],...
              ['C' num2str(ic) ',  \phi_0 = ' num2str(phi0) ' V,  offs_{\perp} = ' num2str(peroffs) ' mV/m,  offs_{||} = ' num2str(paroffs) ' mV/m'],...
              ['E_{||} = E_{||} + ' num2str(R) '*E_{\perp}']};
    title(h(1),strTitle)
end


if 0
    P = cell(1,6)
    Ptmp = zeros(1,nz);    
    for ii=1:6;
        for jj=2:nz
            Ptmp(jj) = trapz(z(1:jj),Ez{1}(1:jj)+ii*Er{1}(1:jj));
        end
        P{ii} = Ptmp;
    end
    
    figure(62)
    set(gcf,'position',[560 322 463 626]);
    h(1) = subplot(2,1,1); h(2) = subplot(2,1,2);
    hca = h(1); 
    plot(z,Ez{1}+1*Er{1},z,Ez{1}+2*Er{1},z,Ez{1}+3*Er{1},z,Ez{1}+4*Er{1},z,Ez{1}+5*Er{1},z,Ez{1}+6*Er{1},z,Er{1}); 
    xlabel(hca,'z [km]'); ylabel(hca,'E_z'); set(hca,'xlim',[-zl zl]); hold(hca,'on');
    hca = h(2);
    
    
    
    
    
end
