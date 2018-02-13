switch ic
  case {11,22} % good model for mms1 and mms2
    limN = 30e3;
    Er = 1*2*1*0.5*1e-3; % reconnection electric field, V/m
    Ey_inflow = 0.5e-3;
    E0 = 1*2*-1e-3;
    Enas = 1*-1.5e-3;
    offsEnas = 15e3;
    B0 = 10e-9; % asymptotical magnetic field, T
    d = 16e3; % thickness of current sheet, m
    dE = 12e3; % thickness of current sheet, m
    b = 8e3;8e3; % bifurcation length scale, m
    dzE = 6e3;
    Bg = 1*4.5*1e-9; % guide field, T
    Bn = 2.0*1e-9; % normal field, T
    BH = 7e-9;
    BLoffset = 0e-9;
    EL = 0*-0.5*1e-3;
    zoffsy = 1e3;
  case {1,2} % good model for mms1 and mms2
    lscale = 1;100/70; % make everything wider
    a = 11e3*lscale; % thickness of current sheet, m
    b = 8e3*lscale; % bifurcation length scale, m
    c = 10e3;
    d = 15e3;
    g = 11e3*lscale;
    
    limN = 30e3;
    Er = -1.0*1e-3; % reconnection electric field, V/m
    Ey_inflow = 0*-0.2e-3;
    E0 = 1*2*-1e-3;
    Enas = -2.5e-3;
    offsEnas = 15e3;
    B0 = 10e-9; % asymptotical magnetic field, T
    d = 8e3; % thickness of current sheet, m
    dE = 12e3; % thickness of current sheet, m
    b = 8e3; % bifurcation length scale, m
    dzE = 6e3;
    Bg = 1*4.5*1e-9; % guide field, T
    Bn = 2.0*1e-9; % normal field, T
    BH = 5e-9;
    BLoffset = 0e-9;
    EL = 0*-0.5*1e-3;
    zoffsy = 1e3;
    %% for vcs = 100 km/s, multiply all lengths with 100/70 = 1.4
    lscale = 1;100/70; % make everything wider
    
    a = 11e3*lscale; % thickness of current sheet, m
    b = 8e3*lscale; % bifurcation length scale, m
    c = 8e3*lscale;
    d = 13e3*lscale;
    g = 12e3*lscale;
    
    zoffsy = 1*1e3*lscale;
    
    
    limN = 30e3*1e-4*lscale;
    Er = -1.0*1e-3; % reconnection electric field, V/m
    Ey_inflow = 0*-1e-3;
    E0 = 1*2*-1e-3;
    Enas = -2.5e-3;
    offsEnas = 15e3*lscale;
    B0 = 10e-9; % asymptotical magnetic field, T
    
    dzEr = 5e3;
    dE = 12e3*lscale; % thickness of current sheet, m
    
    dzE = 6e3*lscale;
    Bg = 1*5*1e-9; % guide field, T
    Bn = 2*1e-9; % normal field, T
    BH = 4e-9;
    BLoffset = 0e-9;
    EL = 0*-0.5*1e-3;            
  case {3,4} % good model for mms2 and mms3
    limN = 30e3;
    Er = 0*-1e-3; % reconnection electric field, V/m
    Ey_inflow = 0*-0.5e-3;
    E0 = -0.5*6e-3;
    B0 = 10e-9; % asymptotical magnetic field, T
    d = 12e3; % thickness of current sheet, m
    dE = 12e3; % thickness of current sheet, m
    b = 6e3; % bifurcation length scale, m
    dzE = 6e3;
    Bg = 5*1e-9; % guide field, T
    Bn = 2.0*1e-9; % normal field, T
    BH = 5e-9;
    BLoffset = 1.5e-9;
    EL = 0*-0.5e-3;    
    zoffsy = 3e3;
    Enas = -0e-3;
    offsEnas = 0;
    a = 11e3*lscale; % thickness of current sheet, m
    b = 8e3*lscale; % bifurcation length scale, m
    c = 10e3;
    d = 15e3;
end

%Bx = @(z) -B0*tanh(z*pi/d);
Bx = @(x,y,z) x*0 + y*0 - 0*1e-9 - 0*abs(z)/a*0.05*B0+- B0*tanh(z/a).*(1-exp(-z.^2/(2*b^2)))-BLoffset; % L
By = @(x,y,z) x*0 + y*0 + 0*z/d*0.03*B0+ Bg-4*BH*sin((z-zoffsy)/c).*exp(-(z-zoffsy).^2/(d^2)).*(1-exp(-(z-zoffsy).^2/(2*b^2))); % M
%By = @(x,y,z) x*0 + y*0 + z*0 + Bg-7*BH*tanh(z/d).*exp(-z.^2/d^2).*(1-exp(-z.^2/(2*b^2))); 
Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn; % N

% new model
if 1
  Bx = @(x,y,z) x*0 + y*0 - 4.5*1e-9*tanh((z+11*1e3)/(5*1e3))- 5.5*1e-9*tanh((z-11*1e3)/(5*1e3)) - 1.5*1e-9; % L
  By = @(x,y,z) x*0 + y*0 + Bg + 5.5*1e-9*exp(-(z+11*1e3).^2/((5*1e3)^2)) - 4.5*1e-9*exp(-(z-14*1e3).^2/((8*1e3)^2)); % M
  Bz = @(x,y,z) x*0 + y*0 + z*0 + 1.5*1e-9; % N
end

Ex = @(x,y,z) x*0 + y*0 + z*0 + 0*EL*(-exp(-(z.^2)/d^2));
Ey = @(x,y,z) x*0 + y*0 + z*0 + 0*(Er*(exp(-((z-dzEr).^2)/g^2))) + Ey_inflow;
%Ez = @(x,y,z) x*0 + y*0 + z*0 + 1*(E0*sin((z-dzE)/dE).*exp(-(z-dzE).^2/(2*b^2)) + 1*Enas*exp(-(z-offsEnas).^2./d.^2) + -0.8*Enas*exp(-(z+offsEnas).^2./2/d.^2));
Ez = @(x,y,z) x*0 + y*0 + z*0 + 0*(2e-3*exp(-(z+0*1e3).^2/((7*1e3)^2)) - 2*1e-3*exp(-(z-10*1e3).^2/((7*1e3)^2))); % N

B = @(x,y,z) sqrt(Bx(x,y,z).^2 + By(x,y,z).^2 + Bz(x,y,z).^2);
%% integrate magnetic field from z = -30 to z = 30 to see how it looks and
if 0
  %%
% so that one can overplot in figures
zBstart = -30*1e3;
zBstop = 30*1e3;
%xyzB(1,:) = [zBstart 0 0];
nz = 5; % Bn is constant, so step in this
xyz = zeros(nz,3);
xyz(:,3) = linspace(zBstart,zBstop,nz);
dz = (zBstop-zBstart)/nz;
dr = dz;
dr = 100;
if 0
for iz = 2:nz  
    dx = dz*Bx(xyz(iz-1,1),xyz(iz-1,2),xyz(iz-1,3))/Bz(xyz(iz-1,1),xyz(iz-1,2),xyz(iz-1,3));
    dy = dz*By(xyz(iz-1,1),xyz(iz-1,2),xyz(iz-1,3))/Bz(xyz(iz-1,1),xyz(iz-1,2),xyz(iz-1,3));
    %fprintf('dx = %g  dy = %g \n',dx,dy)
    xyz(iz,1) = xyz(iz-1,1) + dx;
    xyz(iz,2) = xyz(iz-1,2) + dy;  
end
else
  xyz = xyz(1,:);
  while xyz(end,3)<zBstop
    dx = dr*Bx(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    dy = dr*By(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    dz = dr*Bz(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    %sqrt(dx.^2 + dy.^2 +dz.^2)
    %fprintf('dx = %g  dy = %g   dy = %g \n',dx,dy,dz)
    xyz(end+1,:) = xyz(end,:) + [dx dy dz];    
    %fprintf('xyz(end) = %g %g %g \n',xyz(end,1),xyz(end,2),xyz(end,3))
  end
end


if 0
  %zz = tocolumn(linspace(zBstart,zBstop,nBz)); 
%xyzB = cumsum([Bx(zz*0,zz*0,zz) By(zz*0,zz*0,zz) Bz(zz*0,zz*0,zz)],1);
hca = subplot(1,1,1);
% plot3(hca,xyz(:,1),xyz(:,2),xyz(:,3));
% hca.XLabel.String = 'X';
% hca.YLabel.String = 'Y';
% hca.ZLabel.String = 'Z';

if 0 % LMN
  plot3(hca,xyz(:,3)*1e-3,xyz(:,1)*1e-3,xyz(:,2)*1e-3);
  hca.XLabel.String = 'N'; % Z
  hca.YLabel.String = 'L'; % X
  hca.ZLabel.String = 'M'; % Y
elseif 0 % L-M-N
  plot3(hca,xyz(:,3)*1e-3,-xyz(:,2)*1e-3,xyz(:,1)*1e-3);
  hca.XLabel.String = 'N'; % Z
  hca.YLabel.String = '-M'; % -Y
  hca.ZLabel.String = 'L'; % X
elseif 1 % L-M-N
  
  qscale = 0.5;
  plot3(hca,xyz(:,3)*1e-3,-xyz(:,2)*1e-3,xyz(:,1)*1e-3);
  hold(hca,'on')
  nq = 30;
  toplot = fix(linspace(1,size(xyz,1),nq));
  quiver3(hca,xyz(toplot,3)*1e-3,-xyz(toplot,2)*1e-3,xyz(toplot,1)*1e-3,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
    -By(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale)
  
  hca.XLabel.String = 'N'; % Z
  hca.YLabel.String = '-M'; % -Y
  hca.ZLabel.String = 'L'; % X
else % L-M-N
  xyz = xyz200; 
  qscale = 0.5;
  plot3(hca,xyz(:,3)*1e-3,-xyz(:,2)*1e-3,xyz(:,1)*1e-3);
  hold(hca,'on')
  nq = 30;
  toplot = fix(linspace(1,size(xyz,1),nq));
  quiver3(hca,xyz(toplot,3)*1e-3,-xyz(toplot,2)*1e-3,xyz(toplot,1)*1e-3,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
    -By(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale)
   
  xyz = xyz1000; plot3(hca,xyz(:,3)*1e-3,-xyz(:,2)*1e-3,xyz(:,1)*1e-3);
  toplot = fix(linspace(1,size(xyz,1),nq));
  quiver3(hca,xyz(toplot,3)*1e-3,-xyz(toplot,2)*1e-3,xyz(toplot,1)*1e-3,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
    -By(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale)
  
  xyz = xyz2000; plot3(hca,xyz(:,3)*1e-3,-xyz(:,2)*1e-3,xyz(:,1)*1e-3);
  toplot = fix(linspace(1,size(xyz,1),nq));
  quiver3(hca,xyz(toplot,3)*1e-3,-xyz(toplot,2)*1e-3,xyz(toplot,1)*1e-3,...
     Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
    -By(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale)
  
  hca.XLabel.String = 'N'; % Z
  hca.YLabel.String = '-M'; % -Y
  hca.ZLabel.String = 'L'; % X
end
axis(hca,'equal')

if 0 % add quivers to see if it's ok  
  dxyz_ = 1; % 1 km = 1000 m
  xx_ = linspace(xyz(end,1),xyz(end,1),fix(abs((xyz(1,1)-xyz(end,1))/dxyz_)));
  yy_ = linspace(xyz(end,2),xyz(end,2),fix(abs((xyz(1,2)-xyz(end,2))/dxyz_)));
  zz_ = linspace(xyz(end,3),xyz(end,3),fix(abs((xyz(1,3)-xyz(end,3))/dxyz_))); 
  [XX,YY,ZZ] = meshgrid(xx_,yy_,zz_);
  quiver(hca,XX*1e-3,YY*1e-3,ZZ*1e-3,Bx(XX,YY,ZZ),By(XX,YY,ZZ),By(XX,YY,ZZ))
elseif 0
  xx_ = linspace(xyz(1,1),xyz(end,1),10);
  yy_ = linspace(xyz(1,2),xyz(end,2),10);
  zz_ = linspace(xyz(1,3),xyz(end,3),10); 
  [XX,YY,ZZ] = meshgrid(xx_,yy_,zz_);
  quiver3(hca,ZZ*1e-3,-YY*1e-3,XX*1e-3,Bz(XX,YY,ZZ)*1e9,-By(XX,YY,ZZ)*1e9,Bx(XX,YY,ZZ)*1e9)
elseif 0
  dxyz_ = 3; % 1 km = 1000 m
  xx_ = linspace(xyz(1,1),xyz(end,1),fix(abs((xyz(1,1)-xyz(end,1))/dxyz_)));
  yy_ = linspace(xyz(1,2),xyz(end,2),fix(abs((xyz(1,2)-xyz(end,2))/dxyz_)));
  zz_ = linspace(xyz(1,3),xyz(end,3),fix(abs((xyz(1,3)-xyz(end,3))/dxyz_))); 
  [XX,YY,ZZ] = meshgrid(xx_,yy_,zz_);
  quiver3(hca,ZZ*1e-3,-YY*1e-3,XX*1e-3,Bz(XX,YY,ZZ),-By(XX,YY,ZZ),Bx(XX,YY,ZZ))
end
end


end
%%
if 0 % See what part of the model field that is par and perp
  %%
B_ = [Bx(0,0,zObs*1e3) By(0,0,zObs*1e3) Bz(0,0,zObs*1e3)]; 
Bnorm = irf_norm(B_);
Evec = [Ex(0,0,zObs*1e3) Ey(0,0,zObs*1e3) Ez(0,0,zObs*1e3)];
Epar =  sum(Evec.*Bnorm,2);
Eperp = Evec - Bnorm.*repmat(Epar,1,3);

% model curvature of B
modCurvB = zeros(numel(zObs),3);
dz = zObs(2) - zObs(1);
modCurvB(2:end,:) = repmat(Bnorm(2:end,3),1,3).*diff(Bnorm)/dz; %nT/km

modCurvB_2 = zeros(numel(zObs),3);
dz = zObs(2) - zObs(1);
modCurvB_2(2:end,:) = repmat(B_(2:end,3),1,3).*diff(B_)/dz; %nT/km
end
