f = @(opi,ope,opeb,vti,vte,vteb,om,k) (1 - ...
    opi^2./om^2 - ...
    opeb^2./((om-k*vteb).^2) + ...
    2*ope^2./k.^2/vte^2+ ...
    2i*sqrt(pi)*ope^2*om./k.^3/vte^3);

% Plasma parameters
B=25; n=0.04; no=0;  Te=1600; Ti=2000; Teb = 60; nb=0.02; np = n+nb;
ld = irf_plasma_calc(B,n,no,Te,Ti,'Ld'); % m
ope = irf_plasma_calc(B,n,no,Te,Ti,'Fpe')*2*pi; % rad/s
vte = irf_plasma_calc(B,n,no,Te,Ti,'Vte'); % km/s
vteb = irf_plasma_calc(B,nb,no,Teb,Ti,'Vte'); % km/s
opeb = irf_plasma_calc(B,nb,no,Teb,Ti,'Fpe')*2*pi; % rad/s
opi = irf_plasma_calc(B,np,no,Te,Ti,'Fpp')*2*pi; % rad/s
vti = irf_plasma_calc(B,n,no,Te,Ti,'Vtp'); % km/s


kmin = 0.01/ld;
kmax = 4/ld;
k=logspace(log10(kmin),log10(kmax),100); 
w=[];
x0 = 10+0.1i;
for kk=k
    x=x0;
    options=optimset('Display','off'); 
    x=fsolve(@(om,k)f(ope,ope,opeb,vti,vte,vteb,om,k),x,options,kk); 
    %x=fzero(@(om)f(ope,ope,opeb,vti,vte,vteb,om,kk),x); 
    w=[w,x];
end

wr=real(w);
wi=imag(w);
clear w
plot(k*ld,wr/ope,'o',k*ld,wi/ope,'-')
legend('wr','wi')