n=0.06*1e6; % m^-3
T=1600; % eV
vt=cn_eV2v(T,'eV')*1e3; % m/s
Elim=200; % eV
vlim=cn_eV2v(Elim,'eV')*1e3;  % m/s

f = @(v) exp(-v.^2/vt/vt)/(pi^1.5*vt^3);

vlims = linspace(-vlim,vlim,1000); dvlim = vlims(2)-vlims(1);
vinfs = linspace(-vlim*10,vlim*10,2000); dvinf = vinfs(2)-vinfs(1);
Qlim = trapz(f(vlims))*dvlim;
Q = trapz(f(vinfs))*dvinf;

Qlim^3/(Q^3)




