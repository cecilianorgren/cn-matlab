% Correcting exam
pot1=@(rho,phi)(-2*rho.*cos(phi));
pot2=@(rho,phi,a)(-a^2*cos(phi)./rho-rho.*cos(phi));

a0=0.0; a1=1; a2=5; na2=20; na1=5;
r1=linspace(a0,a1,na1); 
r2=linspace(a1,a2,na2);
nphi=25; phi1=linspace(0.0,2*pi,nphi); phi2=phi1;

[R1,PHI1]=meshgrid(r1,phi1); [R2,PHI2]=meshgrid(r2,phi2);
X1=R1.*cos(PHI1); X2=R2.*cos(PHI2);
Y1=R1.*sin(PHI1); Y2=R2.*sin(PHI2);
Z1=pot1(R1,PHI1); Z2=pot2(R2,PHI2,a1);
surfc(X2,Y2,Z2); hold on; surfc(X1,Y1,Z1); hold off;
ylabel('y');xlabel('x')

%%
[E2x E2y]=gradient(-Z2);
E2r=E2x.*cos(PHI2)+E2y.*sin(PHI2);
E2phi=-E2x.*sin(PHI2)+E2y.*cos(PHI2);
E2y=E2y./R2;

E2rB=E2x.*cos(PHI2)-E2y.*sin(PHI2);
E2phiB=E2x.*sin(PHI2)+E2y.*cos(PHI2);
V=linspace(min(min(Z2)),max(max(Z2)),18);
contour(X2,Y2,Z2,V); hold on;
quiver(X2,Y2,E2rB,E2phiB); 
%
[E1x E1y]=gradient(-Z1); 
E1y=E1y./R1;

E1r=E1x.*cos(PHI1)+E1y.*sin(PHI1);
E1phi=-E1x.*sin(PHI1)+E1y.*cos(PHI1);

E1rB=E1x.*cos(PHI1)-E1y.*sin(PHI1);
E1phiB=E1x.*sin(PHI1)+E1y.*cos(PHI1);

contour(X1,Y1,Z1,V); hold on; 
quiver(X1,Y1,E1rB,E1phiB); hold off;
ylabel('y');xlabel('x')


