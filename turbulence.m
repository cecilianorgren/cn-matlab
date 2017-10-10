t1 = 0; t2 = 999; nt = 100;
t = linspace(t1,t2,nt);

omega1 = 1; omega2 = 100; nomega = 200;
omega = linspace(omega1,omega2,nomega);


L = 100;
x1 = 0; x2 = L; nx = 100; x =linspace(x1,x2,nx);
y1 = 0; y2 = L; ny = 100; y =linspace(y1,y2,ny);

[X,Y] = meshgrid(x,y);

wave = @(omega,t,kx,ky,x,y)(cos(kx*x+ky*y-omega*t));

for ii = 1:nt
    WAVE = wave(omega2,t(ii),2*pi/10,0,X,Y) +...
           wave(omega2,t(ii),0,2*pi/8.2,X,Y);
    surf(X,Y,WAVE); shading flat;
    set(gca,'zlim',[-5 5])
    pause(0.1)
end

%l1 = L;
%k1 =