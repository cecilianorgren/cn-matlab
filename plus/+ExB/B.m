%% Make current matrix
% Integrate the magnetic field everywhere, using cartesian velocities
% Biot-Savart's law
% B = mu0/4pi*intdV(jxr/r^2)
for sim = [0 1]
%sim = 1;
n = 0.07;
mu0 = 1.2566e-6;
e = 1.6022e-19;
rlimproc = 0.6;

% Grid volume
volume = (xGrid(2)-xGrid(1))*(yGrid(2)-yGrid(1))*(zGrid(2)-zGrid(1))*1e9;

% The source velocities
if sim % simulation
    sVAZ = sumMVxyz./sumMxyz; % m/s
    sVX = sumMVXxyz./sumMxyz; % m/s
    sVY = sumMVYxyz./sumMxyz; % m/s
    sVZ = (sumMVZxyz./sumMxyz)*0; % m/s

    sVX = permute(sVX,[2 1 3]);
    sVY = permute(sVY,[2 1 3]);
    sVZ = permute(sVZ,[2 1 3]);
else  % model ExB drift
    [XS,YS,ZS] = meshgrid(xSurf,ySurf,zSurf);
    Ex = XS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
    Ey = YS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
    Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
    sVX = Ey/(B0*1e-9); % m/s 
    sVY = -Ex/(B0*1e-9); % m/s 
    sVZ = sVY * 0;
end

% The source coordinates
[sX,sY,sZ] = meshgrid(xSurf,ySurf,zSurf); % km
sX = sX*1e3; sY = sY*1e3; sZ = sZ*1e3;    % m
sR = sqrt(sX.^2 + sY.^2);                 % m

% The field coordinates
fX = sX; 
fY = sY;
fZ = sZ;

% Put to NaN all values outside a certain radius
rmax = rlimproc*rlim;   % km
rind = find(sR>rmax*1e3);
sVX(rind) = NaN;
sVY(rind) = NaN;
sVZ(rind) = NaN;

% Try to take away some noise
if 0
sV = sqrt(sVX.^2 + sVY.^2 + sVZ.^2);
vhigh = find(sV>2000*1e3);
sVX(vhigh) = NaN;
sVY(vhigh) = NaN;
sVZ(vhigh) = NaN;
end

% The magnetic field matrices
Bx = zeros(size(sX)); 
By = zeros(size(sY));
Bz = zeros(size(sZ));

% The source currents
sJX = -e*n*1e6*sVX; % A/m^2 ?
sJY = -e*n*1e6*sVY;
sJZ = -e*n*1e6*sVZ;

if 1 % Without for-loops. Use VX, VY, VZ
xind = 10:2:numel(xSurf)-10;  
yind = 10:2:numel(ySurf)-10;  
zind = 60;%20:20:numel(zSurf); % fix(numel(zSurf)/2);
xn = 0;
tic
for xx = xind  
    xn = xn + 1; disp(['xx = ' num2str(xn) '/' num2str(numel(xind))])
    for yy = yind        
        if sqrt(sX(xx,yy,1).^2+sY(xx,yy,1).^2)*1e-3<rlim*rlimproc
        for zz = zind                        
            fR = [sX(xx,yy,zz) sY(xx,yy,zz) sZ(xx,yy,zz)];
            X = fX(xx,yy,zz) - sX;
            Y = fY(xx,yy,zz) - sY;
            Z = fZ(xx,yy,zz) - sZ;            
            R = sqrt(X.^2 + Y.^2 + Z.^2);
            if R == 0; continue; end
            % JxR_x = sJY.*Z-sJZ.*Y;
            % JxR_y = sJZ.*X-sJX.*Z;
            % JxR_z = sJX.*Y-sJY.*X;
            newBx = nansum(nansum(nansum((sJY.*Z-sJZ.*Y)./R./R./R*volume*mu0/4/pi)));
            newBy = nansum(nansum(nansum((sJZ.*X-sJX.*Z)./R./R./R*volume*mu0/4/pi)));
            newBz = nansum(nansum(nansum((sJX.*Y-sJY.*X)./R./R./R*volume*mu0/4/pi)));
            Bx(xx,yy,zz) = newBx;            
            By(xx,yy,zz) = newBy;
            Bz(xx,yy,zz) = newBz;            
        end
        end
    end
end
toc
end

if sim; simBx = Bx; simBy = By; simBz = Bz; 
else modBx = Bx; modBy = By; modBz = Bz; end
end

%% Plot V and B
if 0 
% Plot magnetic field in red
quiver3(fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,fZ(xind,yind,zind)*1e-3,...
        Bx(xind,yind,zind)*1,By(xind,yind,zind)*1,Bz(xind,yind,zind),'r'); hold on;
% Plot current fiVeld in blue
quiver3(sX(xind,yind,zind)*1e-3,sY(xind,yind,zind)*1e-3,sZ(xind,yind,zind)*1e-3,...
        sJX(xind,yind,zind),sJY(xind,yind,zind),sJZ(xind,yind,zind),'b'); hold off;
end
%% Compare the magnetic fields
if 1
% Plot magnetic field in red
quiver3(fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,fZ(xind,yind,zind)*1e-3,...
        modBx(xind,yind,zind)*1,modBy(xind,yind,zind)*1,modBz(xind,yind,zind),'r'); hold on;
% Plot magnetic field in blue
quiver3(fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,fZ(xind,yind,zind)*1e-3,...
        simBx(xind,yind,zind)*1,simBy(xind,yind,zind)*1,simBz(xind,yind,zind),'b'); hold off;
end
%% Plot difference
if 1
diffBx = modBx-simBx; diffBy = modBy-simBy; diffBz = modBz-simBz;

quiver3(fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,fZ(xind,yind,zind)*1e-3,...
        diffBx(xind,yind,zind)*1,diffBy(xind,yind,zind)*1,diffBz(xind,yind,zind),'b'); hold off;
end
%% Plot surface plots of B next to each other with same colorbar

fig = figure(13);
for k = 1:6; h(k) = subplot(2,3,k); end
set(fig,'position',[839   515   838   441]); % not full screen
isub = 1;
%irf_colormap('poynting'); 
simB = sqrt(simBx.^2+simBy.^2+simBz.^2);
simBp = sqrt(simBx.^2+simBy.^2);
modB = sqrt(modBx.^2+modBy.^2+modBz.^2);
modBp = sqrt(modBx.^2+modBy.^2);

% Theoretical dB at the origin
tdB = e*phi0*n*1e6*mu0/(B0*1e-9)*art2.g(0.999999*lr/lz)*1e9; % nT

irf_colormap('poynting')
% From simulation
if 1 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,squeeze(fX(xind,yind,zind))*1e-3,squeeze(fY(xind,yind,zind))*1e-3,squeeze(simBz(xind,yind,zind)))
    pcolor(hca,squeeze(fX(xind,yind,zind))*1e-3,squeeze(fY(xind,yind,zind))*1e-3,squeeze(simBz(xind,yind,zind))*1e9);
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'simBz')
    %ch(isub-1)=colorbar('peer',hca);
end
if 1 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(simBp(xind,yind,zind)))
    pcolor(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(simBp(xind,yind,zind))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'simBperp')
end
if 1 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(simB(xind,yind,zind)))
    pcolor(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(simB(xind,yind,zind))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'simB')
end
% From model
if 1 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modBz(xind,yind,zind)))
    pcolor(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modBz(xind,yind,zind))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'modBz')
end
if 1 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modBp(xind,yind,zind)))
    pcolor(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modBp(xind,yind,zind))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'modBperp')
end
if 1 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modB(xind,yind,zind)))
    pcolor(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modB(xind,yind,zind))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'modB')
end

ch=colorbar;
ylabel(ch,'B [nT]')
set(ch, 'Position', [.8514 .11 .0381 .8150])
for ii=1:6
      pos=get(h(ii), 'Position');
      set(h(ii), 'Position', [0.9*pos(1) pos(2) 0.82*pos(3) pos(4)]);
      caxis(h(ii),1*[-1 1]*tdB)
end
ylabel(h(1),'y')
ylabel(h(4),'y')
xlabel(h(4),'x')
xlabel(h(5),'x')
xlabel(h(6),'x')
    
%% Compare the model and simulation velocities

sVAZ = sumMVxyz./sumMxyz; % m/s
sVX = sumMVXxyz./sumMxyz; % m/s
sVY = sumMVYxyz./sumMxyz; % m/s
sVZ = (sumMVZxyz./sumMxyz)*0; % m/s

sVX_sim = permute(sVX,[2 1 3]);
sVY_sim = permute(sVY,[2 1 3]);
sVZ_sim = permute(sVZ,[2 1 3]);

if 1
sV_sim = sqrt(sVX.^2 + sVY.^2 + sVZ.^2);
vhigh = find(sV_sim>3300*1e3);
sVX_sim(vhigh) = NaN;
sVY_sim(vhigh) = NaN;
sVZ_sim(vhigh) = NaN;
end
sV_sim = sqrt(sVX_sim.^2 + sVY_sim.^2);

[XS,YS,ZS] = meshgrid(xSurf,ySurf,zSurf);
Ex = XS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
Ey = YS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
sVX_mod = Ey/(B0*1e-9); % m/s 
sVY_mod = -Ex/(B0*1e-9); % m/s 
sVZ_mod = sVY * 0;
sV_mod = sqrt(sVX_mod.^2 + sVY_mod.^2);
%
% Plot V and V
% Plot magnetic field in red
scale_mod = max(max(max(sVX_mod)))*1e-6;
scale_sim = max(max(max(sVX_sim)))*1e-6;
scale_com = 0;
qscale=1e-5;

%indMax = find(sV_sim==max(max(max(sVX_sim))));
%sVX_mod(indMax) = sVX_sim(indMax);
%sVY_mod(indMax) = sVY_sim(indMax);
%sVZ_mod(indMax) = sVZ_sim(indMax);


% Plot current V field in blue
if 0
h1=quiver3(fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,fZ(xind,yind,zind)*1e-3,...
           sVX_mod(xind,yind,zind),sVY_mod(xind,yind,zind),sVZ_mod(xind,yind,zind)*0,...
           scale_com+scale_mod*0,'b');%,'autoscale','off'); 
hU = get(h1,'UData'); hV = get(h1,'VData'); set(h1,'UData',qscale*hU,'VData',qscale*hV)
end
hold on;
if 1
h2=quiver3(fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,fZ(xind,yind,zind)*1e-3,...
           sVX_sim(xind,yind,zind),sVY_sim(xind,yind,zind),sVZ_sim(xind,yind,zind)*0,...
           scale_com+scale_sim*0,'r');%,'autoscale','off'); 
hU = get(h2,'UData'); hV = get(h2,'VData'); set(h2,'UData',qscale*hU,'VData',qscale*hV)
end



hold off;
view([0 0 1])
legend('model','simulation')
%title('Velocity vectors, sim 10 times smaller in scale..')
set(gca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])