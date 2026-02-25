function vout = vout_sweetparker(B,n)

eps0 = 8.8542e-12;
mu0 = 4*pi*1e-7;
c = 299792458;
mp = 1.6726e-27; % assuming protons contribute to the plasma mass density

m = mp;

rho = m*n*1e6; % kg/m^3


vout = B*1e-9/sqrt(mu0*rho);
%vout = sqrt(B1*1e-9*B2*1e-9*B1+B2)*sqrt(B1+B2)/sqrt(rho1*B2+rho2*B1)/sqrt(mu0);
%vout = sqrt(B1*1e-9*B2*1e-9*B1+B2)*sqrt(B1+B2)/sqrt(rho1*B2+rho2*B1);

disp(['v_out = ' num2str(vout*1e-3,'%.f') ' km/s'])
