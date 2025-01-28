% function [wia,wga,kmax,maxfreq,maxgamma,vph,residual] = dg_solver(vd1,vd2,veth1,veth2,vith,wpe,wpe1,wpe2,wpi,kvec,doposwi)

units = irf_units;

Te = 50; % eV
Ti = 500;
ne = 5; % cc
ni = ne;

veth = sqrt(Te*units.eV*2/units.me); % m/s
vith = sqrt(Ti*units.eV*2/units.mp); % m/s

vde = veth*2.05;
vdi = 0;

wpe = sqrt(units.e^2*ne*1e6/units.eps0/units.me); % s^-1
wpi = sqrt(units.e^2*ni*1e6/units.eps0/units.mp); % s^-1

lamDe = veth/wpe/sqrt(2);
kvec = linspace(0.001,0.8,100)/lamDe;
disp(sprintf('wpe = %s, vde = %s, veth = %s, wpi = %s, vith = %s',wpe,vde,veth,wpi,vith))
%%
doposwi = 0;

[wra,wia,kmax,maxfreq,maxgamma,vph,residual] = dg_solver_2species(vde,vdi,veth,vith,wpe,wpi,kvec,doposwi);
plot(kvec*lamDe,[wra;wia]/wpi); legend('w_r','w_i')
set(gca,'ylim',[-2 5],'ygrid','on')
hold('on')
