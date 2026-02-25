% run art2.E3
doClose = 1;
if doClose; close; end
% Load data
doLoad = 1;
cd /Users/Cecilia/Data/BM/20070831

if doLoad
    load eh_model
    % Electric field
    % Epar? - along B
    % Eper? - perp to B, in spin plane
    % EparAC? - along B
    % EperAC? - perp to B, in spin plane
    
    % Position vector
    % dp1 - perp to B, in spin plane, R3-R4
    % dp2 
    % dp3 - along B  
    
    % Angle between spin plane and B field
    % eb_anglong3, eb_anglong4
end

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
        phi0 = 500;        
        lr = 12;
        lz = 5;
        r = 6.3;
        t0 = toepoch([2007 08 31 10 17 51.794]);
        tint3 = t0 + [-0.03 0.03];
        t0 = toepoch([2007 08 31 10 17 51.728]);
        tint4 = t0 + [-0.04 0.04];
        x3 = 12;
        y3 = -12;
        v = 436;
        R=0.1;

    case 16
        ic = 3;
        peroffs = 2;
        paroffs = 0; 
        phi0 =1250;        
        lr = 18;
        lz = 6;
        r = 6.3;
        t0 = toepoch([2007 08 31 10 17 51.794]);
        tint3 = t0 + [-0.03 0.03];
        t0 = toepoch([2007 08 31 10 17 51.728]);
        tint4 = t0 + [-0.04 0.04];
        x3 = 35;
        y3 = -18;
        v = 436;
        R=-0.0;
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

tint = [min([tint3 tint4]) max([tint3 tint4])];
% Get dx, dy, dz for indicated time interval
dx=cn.mean(irf_tlim(dp1,tint),1); % SP
dy=cn.mean(irf_tlim(dp2,tint),1); % BxSP
dz=cn.mean(irf_tlim(dp3,tint),1); % B

% x3 y3 is defined by model, and 

x4 = x3 + dx;
y4 = y3 + dy;
zl = 3*lz; nz = 100; z3 = linspace(-zl,zl,nz); z4 = z3;
%
[Ex3,Ey3,Ez3] = art2.E3(x3,y3,z3,phi0,lr,lz);
[Ex4,Ey4,Ez4] = art2.E3(x4,y4,z4,phi0,lr,lz);
Er3 = sqrt(Ex3.^2+Ey3.^2);
Er4 = sqrt(Ex4.^2+Ey4.^2);
xylim = 2*max(abs([dx dy]))*[-1 1];
[XS,YS] = cn.meshgrid(linspace(xylim(1),xylim(2),100),linspace(xylim(1),xylim(2),100));
Ex = XS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(0/lz).^2)*1e1; % mV/m
Ey = YS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(0/lz).^2)*1e1; % mV/m       
Er = sqrt(Ex.^2+Ey.^2);
        

% Make figure
set(gcf,'defaultAxesFontSize',16);
set(gcf,'defaultTextFontSize',16);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(gcf,'position',[1000 0 360 950])
for kk = 1:3 h(kk) = subplot(3,1,kk); end
isub = 1;
if 1 % perpendicular view
    hca = h(isub); isub = isub + 1;

    pcolor(hca,XS,YS,Ex); hold(hca,'on'); shading(hca,'flat');
    hh = plot(hca,x3,y3,'go',x4,y4,'bo'); hold(hca,'off');    
    set(hh,'markersize',14,'linewidth',2,'marker','square')%,'ydata',y*1.2*[-1 1],'zdata',z*1.2*[-1 1])            
    set(hca,'xlim',xylim,'ylim',xylim)
    axis(hca,'square');
    xlabel(hca,'\perp 2,sp')
    ylabel(hca,'\perp 1')
    zlabel(hca,'||,B')
    view(hca,[0 0 1])

    text(hca,x3,y3,'      C3')
    text(hca,x4,y4,'      C4')
    hc = colorbar('peer',hca);    
end
if 1 % Cuts of Epar and Eper as eh goes by C3
    hca = h(isub); isub = isub + 1;
    plot(hca,z3,Ez3+R*Er3,'r'); hold(hca,'on');
    plot(hca,z3,Ex3,'r'); 
    xlabel(hca,'z [km]'); 
    ylabel(hca,'E3_z'); 
    ic = 3;
    c_eval('Epar = irf_tlim(Epar?,tint?);',ic);
    c_eval('Eper = irf_tlim(Eper?,tint?);',ic);
    zE = -Epar(:,1)*v; zE = zE - mean(zE);
    plot(hca,zE,Epar(:,2)-paroffs,'g'); 
    plot(hca,zE,Eper(:,2)-peroffs,'g'); 
    zElim3 = sort(zE([1 end]));
    set(hca,'xlim',zElim3);     
    hold(hca,'off')
end
if 1 % Cuts of Epar and Eper as eh goes by C4
    hca = h(isub); isub = isub + 1;
    plot(hca,z4,Ez4+R*Er4,'r'); hold(hca,'on');
    plot(hca,z4,Ex4,'r'); 
    xlabel(hca,'z [km]'); 
    ylabel(hca,'E4_z'); 
    ic = 4;
    c_eval('Epar = irf_tlim(Epar?,tint?);',ic);
    c_eval('Eper = irf_tlim(Eper?,tint?);',ic);
    zE = -Epar(:,1)*v; zE = zE - mean(zE);
    plot(hca,zE,Epar(:,2)-paroffs,'b'); 
    plot(hca,zE,Eper(:,2)-peroffs,'b'); 
    zElim4 = sort(zE([1 end]));
    set(hca,'xlim',zElim4);  
    hold(hca,'off')
end
zElim = [min([zElim3; zElim4]) max([zElim3; zElim4])];
set(h(2),'xlim',zElim); set(h(3),'xlim',zElim);  

    strTitle={['l_{\perp} = ' num2str(lr) ' km,  l_{||} = ' num2str(lz) ' km,  v = ' num2str(v) ' km/s'],...
              ['t_0 = ' irf_time(t0,'isoshort')],...
              ['C' num2str(ic) ',  \phi_0 = ' num2str(phi0) ' V,  offs_{\perp} = ' num2str(peroffs) ' mV/m,  offs_{||} = ' num2str(paroffs) ' mV/m'],...
              ['E_{||} = E_{||} + ' num2str(R) '*E_{\perp}']};
    title(h(1),strTitle)
     
%%    
if 0
 % Add a part of Eperp to Epar, error coming in while constructing Epar from Esp    
    % Figure
      
    
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


%%
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
elseif 0
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
end