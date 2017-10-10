% use disp.hot_buneman
% use script art2.whamp to get indata perhaps

vte1 = irf_plasma_calc(25,0.06,0,t(1)*1e3,2.2,'Vte'); % m/s
vte2 = irf_plasma_calc(25,0.06,0,t(2)*1e3,2.2,'Vte'); % m/s
vti = irf_plasma_calc(25,0.06,0,t(1)*1e3,t(3)*1e3,'Vtp'); % m/s
ope1 = irf_plasma_calc(25,n(1)*1e-6,0,t(1)*1e3,2.2,'Fpe')*2*pi; % rad/s
ope2 = irf_plasma_calc(25,n(2)*1e-6,0,t(2)*1e3,2.2,'Fpe')*2*pi; % rad/s
opi = irf_plasma_calc(25,n(3)*1e-6,0,t(3)*1e3,2.2,'Fpp')*2*pi; % rad/s
vde1 = 0;
vde2 = vte2*1;
vdi = 0;

lamDe = irf_plasma_calc(25,0.06,0,t(1)*1e3,2.2,'Ld'); % m/s background temperature but total density

fr_SI = fr/fn*ope1; % rad/s
fim_SI = fim/fn*ope1; % rad/s
k_SI = z'/zn*lamDe;

ind = 5;
f_in = fr_SI(ind)+1i*fim_SI(ind);
k_in = k_SI(ind);
%%
disp = disp.hot_buneman(f_in,k_in,ope1,ope2,opi,vte1,vte2,vti,vde1,vde2,vdi);
%disp = disp.hot_buneman(fr_SI(ind)+1i*fim_SI(ind),k_SI(ind),ope1,ope2,opi,vte1,vte2,vti,vde1,vde2,vdi);
