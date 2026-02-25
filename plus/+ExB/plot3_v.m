if 1
    [XS,YS,ZS] = meshgrid(xSurf,ySurf,zSurf);
    Ex = XS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
    Ey = YS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
    Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
    sVX = Ey/(B0*1e-9); % m/s 
    sVY = -Ex/(B0*1e-9); % m/s 
    sVZ = sVY * 0;
end
% perpendicular velocity
sVP = sqrt(sVX.^2+sVY.^2); 
% azimuthal velocity
xg0 = find(sX>=0); % right half xy-plane
xs0 = find(sX<0); % left half xy-plane
th = zeros(size(sX));
th(xg0) = atand(sY(xg0)./sX(xg0));
th(xs0) = atand(sY(xs0)./sX(xs0)) + 180;        
sVAZ = -sVX.*sind(th) + sVY.*cosd(th);

%normsVAZ = sVAZ/max(max(max(sVAZ)));

normsVAZ = sVAZ/max(max(max(sVAZ)));
normsVAZ(isnan(normsVAZ))=0;
normsVAZ = sVP/max(max(max(sVP)))/3;
if 1 % Plot arrows    
    pxind = 5:5:numel(xSurf);
    pyind = 5:5:numel(ySurf);
    pzind = 5:5:numel(zSurf);
    quiver3(sX(xind,yind,zind),sY(xind,yind,zind),sZ(xind,yind,zind),...
           sVX(xind,yind,zind),sVY(xind,yind,zind),sVZ(xind,yind,zind))
end
if 0 % plot patches
[spX spY spZ] = sphere(20); spC = spX*0+1; %ones(size(spX));
radius = 5;xSurf(2)-xSurf(1); spX=spX*radius; spY=spY*radius; spZ=spZ*radius;

% Plot vaz
pxind = 5:5:numel(xSurf);
pyind = 5:5:numel(ySurf);
pzind = 0+30:5:numel(zSurf)-30;
tn=0;
for xx = pxind    
    tn=tn+1;
    disp([num2str(tn) '/' num2str(numel(pxind))])
    for yy = pyind
        for zz = pzind
            face = normsVAZ(xx,yy,zz);
            if face > 0.05 && face < 1
            %disp(num2str([sX(xx,yy,zz)*1e-3 sY(xx,yy,zz)*1e-3 sZ(xx,yy,zz)*1e-3 face]))
            hs = surf(spX+sX(xx,yy,zz)*1e-3,spY+sY(xx,yy,zz)*1e-3,spZ+sZ(xx,yy,zz)*1e-3,spC); 
            set(hs,'facealpha',face)
            shading flat; hold(gca,'on');
            colorbar;
            caxis(gca,[-1 2])
            set(gca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1],'zlim',zlim*[-1 1])
            end
            %alpha(hs,normVAZ(xx,yy,zz))
        end
    end
end
end

    