vperp = @(r,pth,wc,vr) (pth./r + wc/2).^2 + vr^2;

vt = 1;
vr = 1;
pth = 1;
wc = 1;
r = linspace(0,2*vt,100);

semilogy(r,    vperp(r,pth,wc,vr),...
         vt/wc,vperp(vt/wc,pth,wc,vr),'*')