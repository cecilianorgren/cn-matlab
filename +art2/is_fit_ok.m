%% make different comparisons between Ex and Ez
Ex = @(x,y,z,phi0,lr,lz) x./lr./lr.*phi0.*exp(-x.^2/2./lr./lr-y.^2/2./lr./lr-z.^2/2./lz./lz);
Ey = @(x,y,z,phi0,lr,lz) y./lr./lr.*phi0.*exp(-x.^2/2./lr./lr-y.^2/2./lr./lr-z.^2/2./lz./lz);
Ez = @(x,y,z,phi0,lr,lz) z./lr./lr.*phi0.*exp(-x.^2/2./lr./lr-y.^2/2./lr./lr-z.^2/2./lz./lz);
%%
phi0 = 500;
lr = 12;
lz = 5;
x=15;
y=10;
z=linspace(-3*lz,3*lz,50);

plot(z,Ey(x,10,z,500,lr,lz),'r',z,Ez(x,10,z,500,lr,lz),'r',...
     z,Ey(x*1.4,1.2*10,z,2*500,lr,lz),'b',z,Ez(x*1.4,1.4*10,z,2*500,lr,lz),'b')




%A = Ex()