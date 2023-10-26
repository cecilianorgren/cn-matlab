syms x y z Pxx(x,y,z) Pxy(x,y,z) Pxz(x,y,z) Pyy(x,y,z) Pyz(x,y,z) Pzz(x,y,z)

r = [x,y,z];
P = [Pxx,Pxy,Pxz;Pxy,Pyy,Pyz;Pxz,Pyz,Pzz];

Pxy = 0;
Pxz = 0;
Pyz = 0;
P1 = [Pxx,Pxy,Pxz];
P2 = [Pxy,Pyy,Pyz];
P3 = [Pxz,Pyz,Pzz];

dP1 = divergence(P1,r);
dP2 = divergence(P2,r);
dP3 = divergence(P3,r);

divP = [dP1,dP2,dP3]
curl(divP,r)

