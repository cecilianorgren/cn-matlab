% use makegrid.m

alpha = pi/2;
phi1 = 0; phi2 = alpha; nphi = 20;
r1 = 0; r2 = 100; nr = 20;

pot = @(r,phi,alpha) r.^(pi/alpha).*sin(pi*phi/alpha);
potoff = @(r) r*0;
Ex = @(r,phi,alpha) -(pi./alpha).*r.^(pi/alpha).*sin(phi.*(pi/alpha-1)) ;
Ey = @(r,phi,alpha) -(pi./alpha).*r.^(pi/alpha).*cos(phi.*(pi/alpha-1)) ;

% 
% EX = Ex(R,PHI,alpha);
% EY = Ey(R,PHI,alpha);
% 
% POT = pot(R,PHI,alpha);
% POTOFF = phioff(R);

%surf(X,Y,POT); shading flat; hold on;
%surf(XOFF,YOFF,PHIOFF); shading flat; hold off;
%quiver(X,Y,EX,EY); hold off;

%
subplot(2,3,1); 
alpha = pi/4;
[R,PHI,X,Y]=makegrid(r1,r2,nr,0,alpha,nphi);
contour(X,Y,pot(R,PHI,alpha)); hold on; 
quiver(X,Y,Ex(R,PHI,alpha),Ey(R,PHI,alpha)); axis equal; hold off;

subplot(2,3,2); 
alpha = 2*pi/4;
[R,PHI,X,Y]=makegrid(r1,r2,nr,0,alpha,nphi);
contour(X,Y,pot(R,PHI,alpha)); hold on; 
quiver(X,Y,Ex(R,PHI,alpha),Ey(R,PHI,alpha)); axis equal; hold off;

subplot(2,3,3); 
alpha = 3*pi/4;
[R,PHI,X,Y]=makegrid(r1,r2,nr,0,alpha,nphi);
contour(X,Y,pot(R,PHI,alpha)); hold on; 
quiver(X,Y,Ex(R,PHI,alpha),Ey(R,PHI,alpha)); axis equal; hold off;

subplot(2,3,4); 
alpha = 4*pi/4;
[R,PHI,X,Y]=makegrid(r1,r2,nr,0,alpha,nphi);
contour(X,Y,pot(R,PHI,alpha)); hold on; 
quiver(X,Y,Ex(R,PHI,alpha),Ey(R,PHI,alpha)); axis equal; hold off;

subplot(2,3,5); 
alpha = 5*pi/4;
[R,PHI,X,Y]=makegrid(r1,r2,nr,0,alpha,nphi);
contour(X,Y,pot(R,PHI,alpha)); hold on; 
quiver(X,Y,Ex(R,PHI,alpha),Ey(R,PHI,alpha)); axis equal; hold off;

subplot(2,3,6); 
alpha = 6*pi/4;
[R,PHI,X,Y]=makegrid(r1,r2,nr,0,alpha,nphi);
contour(X,Y,pot(R,PHI,alpha)); hold on; 
quiver(X,Y,Ex(R,PHI,alpha),Ey(R,PHI,alpha)); axis equal; hold off;



