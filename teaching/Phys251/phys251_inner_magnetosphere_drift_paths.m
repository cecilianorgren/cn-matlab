wE = 2*pi/(24*60*60);
RE = 6371200; % m
me = 9.1094e-31;
mp = 1.6726e-27;
e = 1.6022e-19;
kB = 1.3806e-23;
mu0 = 1.2566e-06;
M = 8.22e22; % A m^2

Er = @(r) -wE*mu0*M./(4*pi*r.^2);



r = RE*linspace(1,15,10);
lon = linspace(0,2*pi,100);

[R,LON] = ndgrid(r,lon);

X = R.*cos(LON);
Y = R.*sin(LON);

ER = Er(R);

EX = ER.*cos(LON);
EY = ER.*sin(LON);


E_conv = 1e7; % V/m
E_conv_kVRE = 0.5;
E_conv = E_conv_kVRE*1e3/RE; %*1e3*RE; % kV/RE - > V/m
EX_conv = 0;
EY_conv = E_conv;


hca = subplot(2,1,1);
scale = NaN;
quiver(hca,X/RE,Y/RE,EX,EY)
hold(hca,'on')
quiver(hca,X/RE,Y/RE,EX*0+EX_conv,EY*0+EY_conv)
quiver(hca,X/RE,Y/RE,EX+EX*0+EX_conv,EY+EY*0+EY_conv)

surf(X/RE,Y/RE,X*0-1,EY+EY_conv)
shading(hca,'flat')
colormap(pic_colors('blue_red'))
%hca.CLim = max(abs(hca.CLim))*[-1 1];
hcb = colorbar(hca);

hold(hca,'off')

hca = subplot(2,1,2);
rvec = RE*linspace(1,15,100);
%plot(hca,rvec/RE,Er(rvec),rvec/RE,rvec*0+E_conv)
semilogy(hca,rvec/RE,abs(Er(rvec)),rvec/RE,rvec*0+abs(E_conv))


%%
% transformation matrix between cartesian and spherical
T = [sin(colat).*cos(lon), sin(colat).*sin(lon), cos(colat);...
     cos(colat).*cos(lon), cos(colat).*sin(lon), -sin(colat);...
     -sin(lon)         , cos(lon)        ,  0];
  
T_sp2cart = T';

Ecart = T_sp2cart*Bsphere;
  
bx = Bcart(1,:);
by = Bcart(2,:);
bz = Bcart(3,:);

%% Potential
e = 1.6022e-19;
kB = 1.3806e-23;
mu0 = 1.2566e-06;
RE = 6.3712e+06; % m
M = 8.22e22; % A m^2
me = 9.1094e-31;
mp = 1.6726e-27;
wE = 2*pi/(24*60*60); % Earth's rotation frequency, rad/s
Econv = 1*[0.5 2]*1e3/(RE); % kV/RE - > V/m (to have V in SI units)


E0 = 0.3*1e-3; % V/m
W_eV = 1000;
W_J = W_eV*e;
%v = sqrt(2*W_J/m);

Bz = @(r) mu0*M./(4*pi*r.^3);
mu = @(r,W_J) W_J./Bz(r);
%q = -e;

phi = @(r,lon,m,q,mu,E0) mu*Bz(r)/q - wE*mu0*M./(8*pi*r) - E0*r.*sin(lon);
phi_gradB = @(r,lon,m,q,mu) mu*Bz(r)/q;
phi_corot = @(r,lon,m,q,W_J) - wE*mu0*M./(8*pi*r);
phi_conv = @(r,lon,m,q,W_J,E0) - E0*r.*sin(lon);


r = RE*linspace(1,15,50);
r = RE*logspace(log10(1),log10(15),120);
lon = linspace(0,2*pi,100);

[R,LON] = ndgrid(r,lon);

X = R.*cos(LON);
Y = R.*sin(LON);



clear h
nrows = 1; ncols = 4; ipanel = 0;
for irow = 1:nrows, for icol = 1:ncols, ipanel = ipanel + 1; h(ipanel) = subplot(nrows,ncols,ipanel); end, end
isub = 1;

% Define particle parameters
Ncont = 30;
R_Bref = 5*RE;
B_ref = Bz(R_Bref);
m = me; q = -e; W_eV = 0e4; mu = W_eV*e/B_ref; E0 = 0.8e-3;
if 1 % phi_gradB
  hca = h(isub); isub = isub + 1;
  PHI = phi_gradB(R,LON,m,q,mu);
  %contourf(hca,X/RE,Y/RE,PHI,Ncont)
  pcolor(hca,X/RE,Y/RE,PHI*1e-3)
  shading(hca,'flat')
  hca.YLabel.String = 'y^{GSM} (RE)';
  hca.XLabel.String = 'x^{GSM} (RE)';
  axis(hca,'equal')
  hca.Title.String = {'Gradient B drift potential'};%,sprintf('m = %.0f m_e, q = %.0f e, W = %g eV, mu = %g J/T', m/me, q/e, W_eV,mu)};
  hcb = colorbar(hca);
  irf_legend(hca,{sprintf('m = %.0f m_e', m/me),...
    sprintf('q = %.0f e', q/e),...
    sprintf('W = %g eV', W_eV),...
    sprintf('R_{ref for mu} = %.0f R_E', R_Bref/RE),...
    sprintf('B(R_{ref for mu}) = %g T', B_ref),...
    sprintf('mu = %g J/T = %g MeV/G',mu,mu*1e-3/e*1e-5),...
    sprintf('E_{conv} = %g V/m ',E0)},[-0.0 1.13],'color','k','fontsize',14);
  irf_legend(hca,{sprintf('E_{conv} = %g V/m ',E0)},[-0.0 1.23],'color','k','fontsize',14);
end 
if 1 % phi_corot
  hca = h(isub); isub = isub + 1;
  PHI = phi_corot(R,LON,m,q,mu);
  %contourf(hca,X/RE,Y/RE,PHI,Ncont)
  pcolor(hca,X/RE,Y/RE,PHI*1e-3)
  shading(hca,'flat')
  hca.YLabel.String = 'y^{GSM} (RE)';
  hca.XLabel.String = 'x^{GSM} (RE)';
  axis(hca,'equal')
  hca.Title.String = {'Corotation potential'};%,sprintf('m = %.0f m_e, q = %.0f e, W = %g eV, mu = %g J/T', m/me, q/e, W_eV,mu)};
  hcb = colorbar(hca);
end
if 1 % phi_conv
  hca = h(isub); isub = isub + 1;
  PHI = phi_conv(R,LON,m,q,mu,E0);
  contourf(hca,X/RE,Y/RE,PHI*1e-3,Ncont)
  shading(hca,'flat')
  hca.YLabel.String = 'y^{GSM} (RE)';
  hca.XLabel.String = 'x^{GSM} (RE)';
  axis(hca,'equal')
  hca.Title.String = {'Convective potential'};%,sprintf('m = %.0f m_e, q = %.0f e, W = %g eV, mu= %g J/T', m/me, q/e, W_eV,mu)};
  hcb = colorbar(hca);
  irf_legend(hca,{sprintf('E_{conv} = %g V/m ',E0)},[0.02 0.98],'color','k','fontsize',14,'fontweight','bold');
end
if 1 % phi_gradB + phi_corot + phi_conv
  hca = h(isub); isub = isub + 1;
  PHI = phi(R,LON,m,q,mu,E0);
  contourf(hca,X/RE,Y/RE,PHI*1e-3,Ncont)
  %pcolor(hca,X/RE,Y/RE,PHI)
  shading(hca,'flat')
  hca.YLabel.String = 'y^{GSM} (RE)';
  hca.XLabel.String = 'x^{GSM} (RE)';
  axis(hca,'equal')
  hca.Title.String = {'Total potential'};%,sprintf('m = %.0f m_e, q = %.0f e, W = %g eV, mu = %g J/T', m/me, q/e, W_eV,mu)};
  hcb = colorbar(hca);
  hcb.YLabel.String = 'Potential (kV)';

  % Find saddle point  
  [nx,ny,nz] = surfnorm(X,Y,PHI);
  [az,el,rho] = cart2sph(nx,ny,nz);   % find azimuth and elevation
  [~,ix] = max(el(:));  
  [~,ind]=max(abs(nz(:)));
  hold(hca,'on')
  contour(X/RE,Y/RE,PHI,PHI(ind)*[1 1],'k:','linewidth',2)
  hold(hca,'off')
  irf_legend(hca,{sprintf('E_{conv} = %g V/m ',E0)},[0.02 0.98],'color','k','fontsize',14,'fontweight','bold');

end
if 0
hca = h(isub); isub = isub + 1;
m = mp; q = e; W_eV = 1e9;
PHI = phi(R,LON,me,-e,W_eV*e);
contourf(hca,X/RE,Y/RE,PHI,Ncont)
shading(hca,'flat')
hca.YLabel.String = 'y^{GSM} (RE)';
hca.XLabel.String = 'x^{GSM} (RE)';
axis(hca,'equal')
hca.Title.String = sprintf('m = %.0f m_e, q = %.0f e, W = %g eV', m/me, q/e, W_eV);
end

%colormap(pic_colors('candy3'))
colormap(irf_colormap('waterfall'))

hlinks = linkprop(h,{'CLim','XLim','YLim'});
hlinks.Targets(1).CLim = 5e1*[-1 1];
hlinks.Targets(1).XLim = 0.99*10*[-1 1];
hlinks.Targets(1).YLim = 0.99*10*[-1 1];

c_eval('axis(h(?),''square'');',1:numel(h))
c_eval('h(?).FontSize = 14;',1:numel(h))
c_eval('h(?).LineWidth = 1;',1:numel(h))

hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
delete(hb(1:3));
compact_panels(0.01,0.01)

hb(4).Position(1) = hca.Position(1) + hca.Position(3)+0.01;
c_eval('h(?).YLabel = [];',2:4)
c_eval('h(?).YTickLabel = [];',2:4)

%% Drift frequencies
c = 299792458;
e = 1.6022e-19;
kB = 1.3806e-23;
mu0 = 1.2566e-06;
RE = 6.3712e+06; % m
M = 8.22e22; % A m^2
me = 9.1094e-31;
mp = 1.6726e-27;
wE = 2*pi/(24*60*60); % Earth's rotation frequency, rad/s
Econv = 1*[0.5 2]*1e3/(RE); % kV/RE - > V/m (to have V in SI units)



E0 = 0.3*1e-3; % V/m
W_eV = 1000;
W_J = W_eV*e;
%v = sqrt(2*W_J/m);

Bz = @(r) mu0*M./(4*pi*r.^3);
mu = @(r,W_J) W_J./Bz(r);

v_rel = @(W,m) c*sqrt(1-1./(W.*e./(m*c^2)+1).^2); 
%v_rel2 = @(W,m) c*sqrt(1-(m*c^2/(W.*e + m*c^2)).^2); 
v_nonrel = @(W,m) sqrt(2*W.*e./m); 
ode
gamma = @(W,m) 1./sqrt(1-v_rel(W,m).^2/c^2);

w_gyro = @(W,m,r) e*Bz(r)./(m.*gamma(W,m));
T_gyro = @(W,m,r) 2*pi./w_gyro(W,m,r);

r_gyro = @(W,m,r) sqrt(2*W*e/m)./w_gyro(W,m,r);

clear h
nrows = 1; ncols = 1; ipanel = 0;
for irow = 1:nrows, for icol = 1:ncols, ipanel = ipanel + 1; h(ipanel) = subplot(nrows,ncols,ipanel); end, end
isub = 1;

Wvec = logspace(0,9,100); % eV
REvec = linspace(1,8,80)*RE;
[WVEC,REVEC] = ndgrid(Wvec,REvec);

if 1 % Relativistic and non-relativistic speeds
  V_rel = v_rel(Wvec,me);
  V_nonrel = v_nonrel(Wvec,me);
  hca = h(isub); isub = isub + 1;
  plot(hca,Wvec,V_rel/c,Wvec,V_nonrel/c,Wvec,Wvec*0+1,'k--')
  %hca.YTick = 10.^(-7:1:10);
  hca.XScale = 'log';
  hca.YScale = 'log';
  hca.XLabel.String = 'Kinetic energy (eV)';
  hca.YLabel.String = 'Speed (c)';
  %hcb = colorbar(hca);
  %hcb.Label.String = 'Relativistic speed, v/c';
  legend(hca,{'Relativistic','Non-relativsistic','Speed of light'},'box','off')
end
%%

if 1 % Larmor frequency
  V_rel = v_rel(WVEC,me);
  hca = h(isub); isub = isub + 1;
  contourf(hca,REVEC/RE,WVEC,V_rel/c,20)
  %hca.YTick = 10.^(-7:1:10);
  hca.YScale = 'log';
  hca.YLabel.String = 'Kinetic energy (eV)';
  hca.XLabel.String = 'L';
  hcb = colorbar(hca);
  hcb.Label.String = 'Relativistic speed, v/c';
end
if 1 % gamma
  Gamma = gamma(WVEC,me);
  hca = h(isub); isub = isub + 1;
  contourf(hca,REVEC/RE,WVEC,Gamma,20)
  %hca.YTick = 10.^(-7:1:10);
  hca.YScale = 'log';
  hca.YLabel.String = 'Kinetic energy (eV)';
  hca.XLabel.String = 'L';
  hcb = colorbar(hca);
  hcb.Label.String = '\gamma';
end
if 1 % Larmor frequency
  F_L = w_gyro(WVEC,me,REVEC)/(2*pi);
  hca = h(isub); isub = isub + 1;
  [C,h] = contourf(hca,REVEC/RE,WVEC,log10(F_L),20);
  clabel(C,h,'LabelSpacing',200,'Color','k','FontWeight','bold')
  %hca.YTick = 10.^(-7:1:10);
  hca.YScale = 'log';
  hca.YLabel.String = 'Kinetic energy (eV)';
  hca.XLabel.String = 'L';
  hcb = colorbar(hca);
  hcb.Label.String = 'Larmor frequency (Hz)';
end
if 1 % Larmor period
  T_L = T_gyro(WVEC,me,REVEC);
  hca = h(isub); isub = isub + 1;
  contourf(hca,REVEC/RE,WVEC,T_L,20)
  %hca.YTick = 10.^(-7:1:10);
  hca.YScale = 'log';
  hca.YLabel.String = 'Kinetic energy (eV)';
  hca.XLabel.String = 'L';
  hcb = colorbar(hca);
  hcb.Label.String = 'Larmor period (s)';
end

%% Grad-B drift
e = 1.6022e-19;
kB = 1.3806e-23;
mu0 = 1.2566e-06;
RE = 6.3712e+06; % m
M = 8.22e22; % A m^2
me = 9.1094e-31;
mp = 1.6726e-27;
wE = 2*pi/(24*60*60); % Earth's rotation frequency, rad/s
Econv = 1*[0.5 2]*1e3/(RE); % kV/RE - > V/m (to have V in SI units)


E0 = 0.3*1e-3; % V/m
W_eV = [1e2 1e3 1e4 1e5];
W_J = W_eV*e;
q = e;

vgradB = @(W,r,q) -W*12*pi.*r.^2/(q*mu0*M);
vcorot = @(r) -wE*r;

rvec = linspace(0,10,100)*RE;
[R,W] = ndgrid(rvec,W_J);

colors = pic_colors('matlab');
clear h
nrows = 1; ncols = 1; ipanel = 0;
for irow = 1:nrows, for icol = 1:ncols, ipanel = ipanel + 1; h(ipanel) = subplot(nrows,ncols,ipanel); end, end
isub = 1;

if 0 % drift speeds
  hca = h(isub); isub = isub + 1;
  plot(hca,rvec/RE,vcorot(rvec)*1e-3,'k--','linewidth',1)
  hold(hca,'on')
  plot(hca,rvec/RE,vgradB(W,R,q)*1e-3,'linewidth',1)
  hold(hca,'off')
  
  hca.LineWidth = 1;
  hca.FontWeight = 'bold';
  hca.XLabel.String = 'L = r_{eq}/R_E';
  hca.YLabel.String = 'Rotation speed (km/s)';
  hca.FontSize = 13;
  for iW = 1:numel(W_eV), leginp(iW).W = W_eV(iW); leginp(iW).q = q/e; end  
  legs = arrayfun(@(x)(sprintf('v_{grad-B}: W = %g eV, q = %ge',x.W,x.q)), leginp, 'UniformOutput', false);
  legs = {'v_{corot}',legs{:}};
  hleg = legend(hca,legs{:},'box','off','location','best');
  hleg.Box = 'off';
  hca.YTick = -1000:20:1000;
end
if 0 % drift frequency
  hca = h(isub); isub = isub + 1;
  t_corot = abs(2*pi*rvec./vcorot(rvec));
  t_gradB = abs(2*pi*R./vgradB(W,R,q));
  
  plot(hca,rvec/RE,1./t_corot,'k--','linewidth',1)
  hold(hca,'on')
  plot(hca,rvec/RE,1./t_gradB,'linewidth',1)
  hold(hca,'off')
  
  hca.LineWidth = 1;
  hca.FontWeight = 'bold';
  hca.XLabel.String = 'L = r_{eq}/R_E';
  hca.YLabel.String = 'Rotation frequency (Hz)';
  hca.FontSize = 13;
  for iW = 1:numel(W_eV), leginp(iW).W = W_eV(iW); leginp(iW).q = q/e; end  
  legs = arrayfun(@(x)(sprintf('v_{grad-B}: W = %g eV, q = %ge',x.W,x.q)), leginp, 'UniformOutput', false);
  legs = {'v_{corot}',legs{:}};
  hleg = legend(hca,legs{:},'box','off','location','best');
  hleg.Box = 'off';
  %hca.YTick = -1000:20:1000;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end
if 1 % drift periods
  hca = h(isub); isub = isub + 1;
  t_corot = abs(2*pi*rvec./vcorot(rvec));
  t_gradB = abs(2*pi*R./vgradB(W,R,q));
  
  plot(hca,rvec/RE,t_corot/60/60,'k--','linewidth',1)
  hold(hca,'on')
  plot(hca,rvec/RE,t_gradB/60/60,'linewidth',1)
  hold(hca,'off')
  
  hca.LineWidth = 1;
  hca.FontWeight = 'bold';
  hca.XLabel.String = 'L = r_{eq}/R_E';
  hca.YLabel.String = 'Rotation period (h)';
  hca.FontSize = 13;
  for iW = 1:numel(W_eV), leginp(iW).W = W_eV(iW); leginp(iW).q = q/e; end  
  legs = arrayfun(@(x)(sprintf('v_{grad-B}: W = %g eV, q = %ge',x.W,x.q)), leginp, 'UniformOutput', false);
  legs = {'v_{corot}',legs{:}};
  hleg = legend(hca,legs{:},'box','off','location','best');
  hleg.Box = 'off';
  %hca.YTick = -1000:20:1000;
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end


