% Check angular momentum of particles
B0 = 20e-9; % nT
% corresponding to Aaz = r*B0/2

ie = 11;
%Stmp = S(3,1,ie);
Stmp = S(ie);
vpar = Stmp.vpar;
vperp = Stmp.vperp;

t = Stmp.t;
x = Stmp.x;
y = Stmp.y;
z = Stmp.z;
r = sqrt(x.^2 + y.^2);
v_r = gradient(r,t); % km/s , m/s?


Az = r*B0/2;

vx = Stmp.vx;
vy = Stmp.vy;
vz = Stmp.vz;

az = atan2(y,x); % azimuthal location in x,y-plane
omega = gradient(az,t); % rad/s 
v_az = r.*omega; % km/s

%v_az = atan2(vy,vx);

%v_az = sqrt(vperp.^2 - v_rad.^2);


vx_rad = vx.*x./r;
vy_rad = vy.*y./r;
v_rad = sqrt(vx_rad.^2 + vy_rad.^2);

ang_mom = v_rad.*r;
ang_mom2 = v_r.*r;

phase_perp = atan2(v_az,v_r);
vperp_check = sqrt(v_az.^2 + v_r.^2);



rmax = max(r);

xlim = rmax*[-1 1];
ylim = rmax*[-1 1];

rr = logspace(1,log10(rmax),10);
  
h = setup_subplots(3,2);
isub = 1;

hca = h(isub); isub = isub + 1;
scatter(hca,x,y,1,vperp)
hca.XLim = xlim;
hca.YLim = ylim;
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'v_{\perp}';
colormap(hca,pic_colors('blue_gray_red'))
if 1 % draw concentric circles for reference
  hold(hca,'on')
  for irr = 1:numel(rr)
    plot(hca,rr(irr)*cosd(0:360),rr(irr)*sind(0:360),'k-')
  end
  hold(hca,'off')
end
axis(hca,'equal')

hca = h(isub); isub = isub + 1;
toplot = v_az;
scatter(hca,x,y,1,toplot)
hca.XLim = xlim;
hca.YLim = ylim;
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'v_{\theta} (azimuthal angle)';
colormap(hca,pic_colors('blue_gray_red'))
hca.CLim = max(abs(hca.CLim))*[-1 1];
try
hca.CLim = prctile(toplot,99)*[-1 1];
end
if 1 % draw concentric circles for reference
  hold(hca,'on')
  for irr = 1:numel(rr)
    plot(hca,rr(irr)*cosd(0:360),rr(irr)*sind(0:360),'k-')
  end
  hold(hca,'off')
end
axis(hca,'equal')

hca = h(isub); isub = isub + 1;
scatter(hca,x,y,1,v_r)
hca.XLim = xlim;
hca.YLim = ylim;
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'v_{r}';
colormap(hca,pic_colors('blue_gray_red'))
hca.CLim = max(abs(hca.CLim))*[-1 1];
if 1 % draw concentric circles for reference
  hold(hca,'on')
  for irr = 1:numel(rr)
    plot(hca,rr(irr)*cosd(0:360),rr(irr)*sind(0:360),'k-')
  end
  hold(hca,'off')
end
axis(hca,'equal')

if 0 % v_r,v_az
hca = h(isub); isub = isub + 1;
plot(hca,v_r,v_az)
hca.XLabel.String = 'v_{r}';
hca.YLabel.String = 'v_{\theta}';
hca.XLim = 8e6*[-1 1];
hca.YLim = 8e6*[-1 1];
end

if 0 % oemga
hca = h(isub); isub = isub + 1;
toplot = omega;
toplot(omega>prctile(omega,99)) = NaN;
plot(hca,z,toplot)
hca.YLabel.String = 'omega';
end

hca = h(isub); isub = isub + 1;
plot(hca,z*1e-3,phase_perp)

hca = h(isub); isub = isub + 1;
plot(hca,phase_perp/pi,(Stmp.Ekpar+Stmp.Ep)/S(ie).phi0,'.')
hca.XLabel.String = '\xi/\pi';
hca.YLabel.String = 'U_{k||}+\phi';
%hold(hca,'on')
%%
hca = h(isub); isub = isub + 1;
plot(hca,z,r)

hca = h(isub); isub = isub + 1;
plot(hca,z,vx,z,vy,z,v_rad)

hca = h(isub); isub = isub + 1;
plot(hca,z,ang_mom,z,ang_mom2)

hca = h(isub); isub = isub + 1;
plot(hca,t,v_rad,t,v_r)