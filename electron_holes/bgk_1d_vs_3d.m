% Check funxctions from Chen et al.
% Given a certain 3D structure of an ESW/EH, how much does the trapped
% density and trapped phase space density matter? My take at the moment is
% that even if the charge pertubation dn = ni-ne may change significantly,
% the trapped density only changes a little, becase the overall variation
% between trapped and passing are so large. I.e., if n0 = 0.04, the free
% density goes down to about half of that, while the trapped density goes
% up to a half... So the total density is still close to n0 = 0.04.
% Therefore, even if dn = changes by a factor of two, the trapped density
% only changes a little, and therefore the phase space density should also
% only change a little.

% Chen assumes the particles are strongly magnetized (based on observations 
% from the auroral region by Ergun et al. 1998). Therefore, the parallel
% motion can sort of be decoupled from the perpendicular motion. And there
% is no trap bouncing in the perpendicular direction. Therefore, a 3D
% structure of an electron hole can be modeled by many 1D electron holes.
% I.e. for each radius, we have one 1D electron hole. The charge
% perturbation is modified from the strictly 1D Poisson's equation.

%% Chen et. al. 2002
% Set up grid to work with
r = linspace(0,3,1001);
z = linspace(-3,3,2001);
rs = linspace(0,1,1); % perpendicular scale size
lz = linspace(0,1,1); % parallel scale size
phi0 = 1; % potential at center of structure

J0 = @(r) besselj(0,r); % Bessel function of zeroth order
absJ0 = @(r) abs(besselj(0,r)); % Bessel function of zeroth order
l00 = fminsearch(@(r) abs(J0(r)),0);

% Parallel potential profile, Gaussian
phi_par = @(z,lz) phi0.*exp(-0.5*(z./lz).^2); 

% Total potential profile
phi = @(r,z) phi_par(z,lz).*J0(l00*r./rs);

% I need to estimate the difference in trapped phase space density.


%% Chen et al. 2002, 1D
Te = 1;
m = 1;
vte = sqrt(2*Te/m);
lz = linspace(0,4.4,1); % parallel scale size
z = linspace(-5,5,99)*lz;

phi0 = 1; % potential at center of structure
v = linspace(-3,3,1001)*sqrt(phi0/Te);
dv = v(2) - v(1);
[Z,V] = meshgrid(z,v); 

% Parallel potential profile, Gaussian
phi = @(z,lz) (phi0/Te).*exp(-0.5*(z./lz).^2); 

% v in terms of vt, phi in terms of phi0
w = @(v,z,lz) v.^2 - phi(z,lz);

fp = @(v,z,lz) (2/sqrt(pi))*exp(-w(v,z,lz));
%ft = @(phi0,v,z,lz) (2*sqrt(w(v,z,lz))./pi/lz.^2);%.*(1-2*log(-w(v,z,lz)/phi0))...
%            + 2*exp(w(v,z,lz))/sqrt(pi).*(1-erf(sqrt(-w(v,z,lz))));

ft1 = @(phi0,lz,w) (4*sqrt(-w)./pi./lz.^2).*(1-2*log(-4*w/phi0));
ft2 = @(phi0,lz,w) 2*exp(-w)./sqrt(pi).*(1-erf(sqrt(-w)));

%ft1 = @(phi0,v,z,lz) (4*sqrt(w(v,z,lz))./pi/lz.^2);
%ft2 = @(phi0,v,z,lz) 2*exp(w(v,z,lz))/sqrt(pi).*(1-erf(sqrt(-w(v,z,lz))));


W = w(V,Z,lz);
WT = W; WT(W>0) = nan; % erf will give error if W is positive
WP = W; WP(W<0) = nan;

FT1 = ft1(phi0,lz,WT);
FT2 = ft2(phi0,lz,WT);

FT = FT1 + FT2;
FT(W>0) = 0;

FP = fp(V,Z,lz);
FP(W<0) = 0;

F = FP + FT;
% Densities
np = sum(FP,1)*dv;
nt = sum(FT,1)*dv;
n = sum(F,1)*dv;

          
nrows = 3;
ncols = 2;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % phi
  hca = h(isub); isub = isub + 1;
  plot(hca,z,phi(z,lz))
  hca.XLabel.String = 'z';
end
if 1 % w
  hca = h(isub); isub = isub + 1;
  levels = floor(min(W(:))):1:ceil(max(W(:)));
  levels = 10;
  contourf(hca,Z,V,W,levels)
  shading(hca,'flat')
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'v';
end
if 1 % fp + ft
  hca = h(isub); isub = isub + 1;
  pcolor(hca,Z,V,F)
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f_p+f_t';
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'v';
end
if 0 % ft
  hca = h(isub); isub = isub + 1;
  pcolor(hca,Z,V,FT)
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f_t';
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'v';
end
if 0 % fp
  hca = h(isub); isub = isub + 1;
  pcolor(hca,Z,V,FP)
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f_p';
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'v';
end
if 1 % cut through center, to see discontinuitites
  hca = h(isub); isub = isub + 1;
  plot(hca,v,F(:,fix(size(F,2)/2)))
  hca.XLabel.String = 'v';
end
if 1 % cut through center, to see discontinuitites
  hca = h(isub); isub = isub + 1;
  plot(hca,z,F(fix(size(F,1)/2),:))
  hca.XLabel.String = 'z';
end
if 1 % densities
  hca = h(isub); isub = isub + 1;
  plot(hca,z,np,z,nt,z,n)
  hca.XLabel.String = 'z';
  legend(hca,{'passing','trapped','total'})
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end

%% Chen et al. 2002, 3D
Te = 1;
m = 1;
vte = sqrt(2*Te/m);
phi0 = 1; % potential at center of structure
% v is normalized to vt
% phi is normalized to Te

% Set up grid to work with
r = linspace(0,0.5,1);
z = linspace(-3,3,101);
rs = linspace(0,1,1); % perpendicular scale size
lz = linspace(0,1,1); % parallel scale size

% Velocity space
v = linspace(-3,3,1001)*sqrt(phi0/Te);
dv = v(2) - v(1);
[Z,V] = meshgrid(z,v); 

% Potential
J0 = @(r) besselj(0,r); % Bessel function of zeroth order
absJ0 = @(r) abs(besselj(0,r)); % Bessel function of zeroth order
l00 = fminsearch(@(r) abs(J0(r)),0);
kperp = l00/rs;

% Parallel potential profile, Gaussian
phi_par = @(z,lz) phi0.*exp(-0.5*(z./lz).^2); 

% Total potential profile
phi = @(r,z) phi_par(z,lz).*J0(l00*r./rs);

% Energy
w = @(r,z,v) v.^2 - phi(r,z);

% Passing electron phase space density
fp = @(r,z,v) (2/sqrt(pi))*exp(-w(r,z,v));

% Trapped electron phase space density


if 0
% Trapped phase space density
ft1 = @(phi0,lz,w) (4*sqrt(-w)./pi./lz.^2).*(1-2*log(-4*w/phi0));
ft2 = @(phi0,lz,w) 2*exp(-w)./sqrt(pi).*(1-erf(sqrt(-w)));


W = w(V,Z,lz);
WT = W; WT(W>0) = nan; % erf will give error if W is positive
WP = W; WP(W<0) = nan;

FT1 = ft1(phi0,lz,WT);
FT2 = ft2(phi0,lz,WT);

FT = FT1 + FT2;
FT(W>0) = 0;

FP = fp(V,Z,lz);
FP(W<0) = 0;

F = FP + FT;
% Densities
np = sum(FP,1)*dv;
nt = sum(FT,1)*dv;
n = sum(F,1)*dv;
end
  
r1 = 0.0*rs;
r2 = 0.5*rs;
r3 = 0.7*rs;
  
nrows = 3;
ncols = 2;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % phi
  hca = h(isub); isub = isub + 1;  
  plot(hca,z,phi(r1,z),z,phi(r2,z),z,phi(r3,z))
  hca.XLabel.String = 'z';
  legend(hca,{sprintf('r=%.1fr_s',r1/rs),sprintf('r=%.1fr_s',r2/rs),sprintf('r=%.1fr_s',r3/rs)},'box','off')
end
if 1 % w
  hca = h(isub); isub = isub + 1;
  levels = floor(min(W(:))):1:ceil(max(W(:)));
  levels = 10;
  contourf(hca,Z,V,W,levels)
  shading(hca,'flat')
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'v';
end
if 1 % fp + ft
  hca = h(isub); isub = isub + 1;
  pcolor(hca,Z,V,F)
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f_p+f_t';
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'v';
end
if 0 % ft
  hca = h(isub); isub = isub + 1;
  pcolor(hca,Z,V,FT)
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f_t';
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'v';
end
if 0 % fp
  hca = h(isub); isub = isub + 1;
  pcolor(hca,Z,V,FP)
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f_p';
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'v';
end
if 1 % cut through center, to see discontinuitites
  hca = h(isub); isub = isub + 1;
  plot(hca,v,F(:,fix(size(F,2)/2)))
  hca.XLabel.String = 'v';
end
if 1 % cut through center, to see discontinuitites
  hca = h(isub); isub = isub + 1;
  plot(hca,z,F(fix(size(F,1)/2),:))
  hca.XLabel.String = 'z';
end
if 1 % densities
  hca = h(isub); isub = isub + 1;
  plot(hca,z,np,z,nt,z,n)
  hca.XLabel.String = 'z';
  legend(hca,{'passing','trapped','total'})
  hca.XGrid = 'on';
  hca.YGrid = 'on';
end



%% Electron test particles
% To see how the electrons move through the structure, I pass them through
% a defined potential.
% Define potential as a function
% To calculate the electric field from the potential, use symbolic
% expressions
% phi0 = 1
% x/lx -> x, y/ly -> y, z/lz -> z 
units = irf_units;
m = units.me;
q = -units.e;
T = 100;
vt = sqrt(T*units.eV*2/units.me)*1e-3;
Ek = 10;
vk = sqrt(Ek*units.eV*2/units.me)*1e-3;
lr = 2*4e3;
lx = lr/sqrt(2);
ly = lr/sqrt(2);
lz = 4e3;

% From data fit 
lx = 14*1e3;
ly = 11*1e3;
lz = 3*1e3;

lx = 7.4e3;
ly = 12.6e3;
lz = 4.3e3;


phi0 = 237;
B0 = 20e-9;
syms x y z phi rho E
R = [x,y,z];

r = sqrt(x^2 + y^2); % r = r/rs
J0 = besselj(0,r); % Bessel function of zeroth order
absJ0 = abs(besselj(0,r)); % Bessel function of zeroth order
mf_absJ0 = matlabFunction(absJ0);
l00 = fminsearch(@(x) abs(mf_absJ0(x,0)),0);
J0 = besselj(0,l00*r/lr); % redefine bessel function with 'normalized' r

% Parallel potential profile, Gaussian
phi_par = phi0.*exp(-0.5*(z./lz).^2); 

% Total potential profile
phi = phi_par*J0; % Chen 2002
phi = phi_par*exp(-0.5*(x./lx).^2)*exp(-0.5*(y./ly).^2);  % Gaussian
mf_phi = matlabFunction(phi);


% Electric field
E = -gradient(phi,R);

rho = divergence(E,R)*units.mu0;
mf_rho = matlabFunction(rho,'vars',R);


mf_Ex = matlabFunction(E(1));
mf_Ey = matlabFunction(E(2));
mf_Ez = matlabFunction(E(3));

% Magnetic field
A = [0, B0*x, 0];
B = curl(A,R);

mf_Ay = matlabFunction(A(2),'vars',R); 
mf_Bx = matlabFunction(B(1),'vars',R);
mf_By = matlabFunction(B(2),'vars',R);
mf_Bz = matlabFunction(B(3),'vars',R);

xmax = lr*2;
ymax = xmax;
zmax = 2*xmax;
xvec = linspace(-xmax,xmax,11);
yvec = linspace(-ymax,ymax,11);
zvec = linspace(-zmax,zmax,21);
[X,Y,Z] = meshgrid(xvec,yvec,zvec);
xvec2 = linspace(-xmax,xmax,101);

%options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,obj.xi(([1 end]))));
%options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,[190 210]));

%options = odeset('RelTol',1e-6);
inpEx = @(x,y,z) mf_Ex(x,y,z);
inpEy = @(x,y,z) mf_Ey(x,y,z);
inpEz = @(x,y,z) mf_Ez(x,y,z);
inpBx = @(x,y,z) mf_Bx(x,y,z);
inpBy = @(x,y,z) mf_By(x,y,z);
inpBz = @(x,y,z) mf_Bz(x,y,z);
EoM = @(t,xyz) eom(t,xyz,inpEx,inpEy,inpEz,inpBx,inpBy,inpBz,m,q); 
boxedge = [-50 50]*1e3; % terminate integration when electron exits this region in z.
options = odeset('AbsTol',1e-7,'AbsTol',1e-9,'Events',@(tt,xx) myEventBoxEdge(tt,xx,boxedge));

% Initial conditions
x_init = [1,0.5*lr*1e-3,40,0,vt,-vk]*1e3; % m, m/s
T = [0 1];

[t,x_sol] = ode45(EoM,T,x_init,options);
r_sol = sqrt(x_sol(:,1).^2 + x_sol(:,2).^2);
[PKS,LOCS]= findpeaks(r_sol); 
mean_peak_distance = diff(LOCS);
r_center = movmean(r_sol,fix(mean(mean_peak_distance)));

% Calculate total energy U = Ek + Ep
Ep = -units.e*mf_phi(x_sol(:,1),x_sol(:,2),x_sol(:,3));
Ek = 0.5*units.me*(x_sol(:,4).^2 + x_sol(:,5).^2 + x_sol(:,6).^2);
Eky = 0.5*units.me*(x_sol(:,5).^2);
Ekpar = 0.5*units.me*(x_sol(:,6).^2);

% Extract vector potential
Ay = mf_Ay(x_sol(:,1),x_sol(:,2),x_sol(:,3));


h = setup_subplots(3,2);
isub = 1;
if 1
  hca = h(isub); isub = isub + 1;
  plot3(hca,x_sol(:,1)*1e-3,x_sol(:,2)*1e-3,x_sol(:,3)*1e-3)
  % hca.XLim = [-10 10];
  % hca.YLim = [-10 10];
  % hca.ZLim = boxedge*1e-3;
  % axis(hca,'equal')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
  hca.XLabel.String = 'x (km)';
  hca.YLabel.String = 'y (km)';
  hca.ZLabel.String = 'z (km)';
  hold(hca,'on')
  quiver3(hca,X*1e-3,Y*1e-3,Z*1e-3,mf_Ex(X,Y,Z),mf_Ey(X,Y,Z),mf_Ez(X,Y,Z))
  hold(hca,'off')
end
if 1 % r
  hca = h(isub); isub = isub + 1;
  plot(hca,t,r_sol*1e-3,t,r_center*1e-3)
  hca.YLim(1) = 0;
  hca.YLabel.String = 'r (km)';
end
if 1 % U = Ek + Ep
  hca = h(isub); isub = isub + 1;
  plot(hca,t,Ek/units.eV,t,Ep/units.eV,t,(Ek+Ep)/units.eV,t,-units.e*Ay)
  legend(hca,{'E_k','E_p','E_k + E_p','E_{ky}+qA_y'})  
end
if 1 % U = Ek + Ep
  hca = h(isub); isub = isub + 1;
  plot(hca,t,m*x_sol(:,5)+q*Ay)
  legend(hca,{'mv_y + qA_y'})  
end
if 1 % phi and rho
  hca = h(isub); isub = isub + 1;
  plotyy(hca,xvec2,mf_rho(xvec2,0,0),xvec2,mf_phi(xvec2,0,0))
  legend(hca,{'rho','phi'})  
end

function [value, isterminal, direction] = myEventBoxEdge(t, x_vect,boxedge)
  % integration is terminated when value changes sign
  % for this setup, value is initially negative
  value      = (x_vect(3)-boxedge(1))*(x_vect(3)-boxedge(2));        
  isterminal = 1;   % Stop the integration
  direction  = 0;
end
function  x_res = eom(t,x_vect,inpEx,inpEy,inpEz,inpBx,inpBy,inpBz,m,q)
  x_ = x_vect(1);
  y_ = x_vect(2);
  z_ = x_vect(3);
  vx_ = x_vect(4);
  vy_ = x_vect(5);
  vz_ = x_vect(6);      

  Ex = inpEx(x_,y_,z_);
  Ey = inpEy(x_,y_,z_);
  Ez = inpEz(x_,y_,z_);   
  Bx = inpBx(x_,y_,z_);
  By = inpBy(x_,y_,z_);
  Bz = inpBz(x_,y_,z_); 
  
%   it = size(x_sol_all,1);
%   x_sol_all(it+1,1) = Ex;
%   x_sol_all(it+1,2) = Ey;
%   x_sol_all(it+1,3) = Ez;
%   x_sol_all(it+1,4) = 0;
%   x_sol_all(it+1,5) = 0;
%   x_sol_all(it+1,6) = Bz;
%   x_sol_all(it+1,7) = t;

  % Equations to be solved
  x_res = zeros(6,1);
  x_res(1) = vx_; % dx/dt = vx;
  x_res(2) = vy_; % dy/dt = vy;
  x_res(3) = vz_; % dz/dt = vz;
  x_res(4) = (q/m)*(Ex + vy_*Bz - vz_*By);
  x_res(5) = (q/m)*(Ey + vz_*Bx - vx_*Bz);
  x_res(6) = (q/m)*(Ez + vx_*By - vy_*Bx);                                              

end      