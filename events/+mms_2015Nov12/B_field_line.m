ic = 1;

%% Observed data
load(['/Users/cno062/Research/Events/2015-11-12_071854/' 'gseEB.mat'])
%tintload = irf.tint('2015-11-12T07:18:54.00Z/2015-11-12T07:19:45.00Z');
%c_eval('gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tintload);',1:4);
%c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tintload);',1:4);

tint_bcs = irf.tint('2015-11-12T07:19:20.116Z/2015-11-12T07:19:22.136Z');
c_eval('[out?,l?,v?] = irf_minvar(gseB?.tlim(tint_bcs));',1:4)
c_eval('mva_all(:,:,?) = v?;',1:4)
mva_mean = v1;mean(mva_all,3);
L = mva_mean(1,:); M = mva_mean(2,:); N = mva_mean(3,:);
coordLabels = {'L','M','N'};
lmn = [L;M;N];

c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';')
c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';')

units = irf_units;

toffset1 = 0; toffset2 = -0.03; toffset3 = 0.1; toffset4 = 0.1;
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
c_eval('tintObs = tintObs+-toffset?;',ic) % tinObs(2) correspond
CS_normal_velocity = 70; % km/s
tintObs = tintObs + 1*[-1 1];

zf0 = 0; % starting/ending point of where to get the f
time = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0/CS_normal_velocity;
c_eval('gseBref = mean(gseB?.tlim(time+[-0.005 0.005]).data,1);',ic)
c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs); obsE = obsE.resample(obsB);'...
],ic)
zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;


%% Integrate magnetic field from model to get magnetic field line
ic = 1;
mms_2015Nov12.Bmodel

zBstart = -30*1e3;
zBstop = 30*1e3;
dr = 100;
xyz = [0 0 zBstart];
while xyz(end,3)<zBstop
  dx = dr*Bx(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  dy = dr*By(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  dz = dr*Bz(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
  xyz(end+1,:) = xyz(end,:) + [dx dy dz];      
end

xyzB = [Bx(xyz(:,1),xyz(:,2),xyz(:,3)) By(xyz(:,1),xyz(:,2),xyz(:,3)) Bz(xyz(:,1),xyz(:,2),xyz(:,3))];

% Obtain electric field along the field line
xyzE = [Ex(xyz(:,1),xyz(:,2),xyz(:,3)) Ey(xyz(:,1),xyz(:,2),xyz(:,3)) Ez(xyz(:,1),xyz(:,2),xyz(:,3))];
xEint = -cumtrapz(xyz(:,1),Ex(xyz(:,1),xyz(:,2),xyz(:,3)));
yEint = -cumtrapz(xyz(:,2),Ey(xyz(:,1),xyz(:,2),xyz(:,3)));
zEint = -cumtrapz(xyz(:,3),Ez(xyz(:,1),xyz(:,2),xyz(:,3)));
xyzEint = [xEint, yEint, zEint];

% Plot 
nrows = 4;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

colors = mms_colors('matlab');
if 1 % xyz vs z
  hca = h(isub); isub = isub + 1;  
  plot(hca,xyz(:,3)*1e-3,[xyz]*1e-3) % sqrt(sum(xyz.^2,2))
  hca.YLabel.String = 'R_{@B} (km)';
  hca.XLabel.String = 'N (km)';
  legend(hca,'L','M','N')
  hca.XLim = [min(xyz(:,3)) max(xyz(:,3))]*1e-3;
end
if 1 % xyzB vs z
  hca = h(isub); isub = isub + 1;
  set(hca,'ColorOrder',colors(1:3,:))
  plot(hca,xyz(:,3)*1e-3,xyzB*1e9,'--')
  hold(hca,'on')
  set(hca,'ColorOrder',colors(1:3,:))  
  plot(hca,zObs,obsB.data)
  hold(hca,'off')
  hca.YLabel.String = 'B_{@B} (nT)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [min(xyz(:,3)) max(xyz(:,3))]*1e-3;
end
if 1 % xyzE vs z
  hca = h(isub); isub = isub + 1;
  set(hca,'ColorOrder',colors(1:3,:))
  plot(hca,xyz(:,3)*1e-3,xyzE*1e3,'--')
  hold(hca,'on')
  set(hca,'ColorOrder',colors(1:3,:))  
  plot(hca,zObs,obsE.data)
  hold(hca,'off')
  hca.YLabel.String = 'E_{@B} (mV/m)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [min(xyz(:,3)) max(xyz(:,3))]*1e-3;
end
if 1 % xyzEint vs z
  hca = h(isub); isub = isub + 1;
  plot(hca,xyz(:,3)*1e-3,xyzEint)
  hca.YLabel.String = '\phi_{@B} (V)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [min(xyz(:,3)) max(xyz(:,3))]*1e-3;
end
if 0 % xyzE vs z
  hca = h(isub); isub = isub + 1;
end

%%

if 0
  %%
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

if 1 % add quivers to see if it's ok  
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

