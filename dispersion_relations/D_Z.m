% ?Surfc ?plot ?of ?Re(z) ?and ?Im(z) ? 
xmin=-2; 
xmax=2; 
dx=0.1;
ymin=-2; 
ymax=2; 
dy=0.1; 
[X,Y]=meshgrid(xmin:dx:xmax,ymin:dy:ymax); 
Zz=faddeeva(complex(X,Y))*1i*sqrt(pi);
surfc(xmin:dx:xmax,ymin:dy:ymax,real(Zz));
axis square;
caxis([-1 1]); 
xlabel('x'); ylabel('y');
title('Re(z)');
figure;surfc(xmin:dx:xmax,ymin:dy:ymax,imag(Zz));
axis square; 
caxis([-1 1]); 
xlabel('x');ylabel('y');
title('Im(z)');

