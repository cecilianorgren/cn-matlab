% Units
units = irf_units;
e = 1.6022e-19;
kB = 1.38e-23;
me = 9.1094e-31;
me = 9.1094e-31;
mp = 1.6726e-27;

% Input parameters
m = me;
n = 1*1e6; % m^(-3)
Tx = 1000;
Ty = 2000;
Tz = 1500;
vdx = 1000e3; % m/s
vdx = 0*0.5*fun_vt(Tx);
vdy = 0;
vdz = 0;

% Define functions
fun_vt = @(T) sqrt(2*units.eV*T/(m)); % m/s
f1D = @(n,v,T,vd) n/((pi)^(1/2)*sqrt(2*units.eV*T/(m)))*exp(-(v-vd).^2./(2*units.eV*T/(m)));
f2D = @(n,vx,Tx,vdx,vy,Ty,vdy) (1/n)*f1D(n,vx,Tx,vdx).*f1D(n,vy,Ty,vdy);
f3D = @(n,vx,Tx,vdx,vy,Ty,vdy,vz,Tz,vdz) (1/n^2)*f1D(n,vx,Tx,vdx).*f1D(n,vy,Ty,vdy).*f1D(n,vz,Tz,vdz);


% Prepare grid
v = fun_vt(max([Tx,Ty,Tz]))*linspace(-3,3,300);
dv = v(2) - v(1);
[VX,VY,VZ] = ndgrid(v,v,v);

% Construct distributions
fxyz = f3D(n,VX,Tx,vdx,VY,Ty,vdy,VZ,Tz,vdz);
fxy_red_z = sum(fxyz,3)*dv;
fx_red_yz = sum(sum(fxyz,3),2)*dv*dv;
mom.n = sum(fxyz,'all')*dv*dv*dv;


% Plot
h = setup_subplots(3,1);
isub = 1; 

hca = h(isub); isub = isub + 1;
plot(hca,v,f1D(n,v,Tx,vdx),v,f1D(n,v,Ty,vdy),v,f1D(n,v,Tz,vdz))
hca.XLabel.String = 'v (m/s)';
legend(hca,{'f(v_x)','f(v_y)','f(v_z)'})

hca = h(isub); isub = isub + 1;
plot(hca,v,f3D(n,v,Tx,vdx,0,Ty,vdy,0,Tz,vdz))
hca.XLabel.String = 'v_x (m/s)';
legend(hca,{'f(v_x,v_y=0,v_z=0)'})

hca = h(isub); isub = isub + 1;
plot(hca,v,fx_red_yz)
hca.XLabel.String = 'v_x (m/s)';
legend(hca,{'\int vy vz f(v_x,v_y,v_z)'})

%% Symbolic expressions
clear all
syms v T vd m n
units = irf_units;

T = 1000*units.eV;
m = units.me;

f = n/(pi)^0.5/sqrt(2*T/m)*exp(-(v-vd)^2/(2*T/m));
fv = symfun(f,v);