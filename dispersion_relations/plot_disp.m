%% several stream instability with thermal ion and electron backgrounds
% densities [cc]
ne1=0.0025;
ne2=0.075;
ne3=0.0005;
ne4=0.0027;
ne5=0.00017;
ni=[ne1+ne2+ne3+ne4+ne5];

n=[ni ne1 ne2 ne3 ne4 ne5];

% plasma frequencies [Hz]
fp=[0.033 9e3 9e3 9e3 9e3 9e3].*sqrt(n);

% temperatures [eV]
T=[3000 15 1700 150 100 10]; 
m=[Units.mp Units.me Units.me Units.me Units.me Units.me];
vt=sqrt(T*Units.eV*2./m)/1000;

% drift velocities [km/s]
vdvt=[0 -0.5 0 2.5 -1.8 -3.5];
vd=vdvt.*vt;

%% dispersion relation
D=@(omega,k)(1-(fp(1)/omega)^2-(fp(2)/(omega-k*vd(2)))^2 ...
        -fpe2^2/(om-k*vde2));
