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


%% Electron test particles
% To see how the electrons move through the structure, I pass them through
% a defined potential.
% Define potential as a function
% To calculate the electric field from the potential, use symbolic
% expressions
% phi0 = 1
% x/lx -> x, y/ly -> y, z/lz -> z 
lr = 10e3;
lz = 10e3;
phi0 = 300;
B0 = 10e-9;
syms x y z
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
phi = phi_par*J0;

% Electric field
E = -gradient(phi,R);

mf_Ex = matlabFunction(E(1));
mf_Ey = matlabFunction(E(2));
mf_Ez = matlabFunction(E(3));

%%
%options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,obj.xi(([1 end]))));
%options = odeset('Events',@exitBox,varargin{:},'Events', @(tt,xx) myEventBoxEdge(tt,xx,[190 210]));

%options = odeset('RelTol',1e-6);
inpEx = @(x,y,z) mf_Ex(x,y,z);
inpEy = @(x,y,z) mf_Ey(x,y,z);
inpEz = @(x,y,z) mf_Ez(x,y,z);
inpBx = @(x,y,z) 0;
inpBy = @(x,y,z) 0;
inpBz = @(x,y,z) 0*B0;
units = irf_units;
m = units.me;
q = -units.e;
EoM = @(t,xyz) eom(t,xyz,inpEx,inpEy,inpEz,inpBx,inpBy,inpBz,m,q); 
boxedge = [-50 50]*1e3; % terminate integration when electron exits this region in z.
options = odeset('AbsTol',1e-7,'AbsTol',1e-9,'Events',@(tt,xx) myEventBoxEdge(tt,xx,boxedge));

% Initial conditions
x_init = [5,20,40,0,0,10*1e3]; % m, m/s
T = [0 0.1];

[t,x_final] = ode45(EoM,T,x_init); % 

h = setup_subplots(1,1);
isub = 1;
hca = h(isub); isub = isub + 1;
plot3(hca,x_final(:,1)*1e-3,x_final(:,2)*1e-3,x_final(:,3)*1e-3)
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
      