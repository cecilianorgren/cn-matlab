x=linspace(0,10,100);
y=x;

[X,Y]=meshgrid(x,y);

stream = @(x,y)(2*x.*y);
velo = @(x,y)(x.^2-y.^2);

subplot(1,3,1);contour(X,Y,stream(X,Y)); title('stream function = 2xy')
subplot(1,3,2);contour(X,Y,velo(X,Y)); title('velocity potential = x^2-y^2')
subplot(1,3,3);contour(X,Y,stream(X,Y)); hold on; contour(X,Y,velo(X,Y)); hold off; title('both')


