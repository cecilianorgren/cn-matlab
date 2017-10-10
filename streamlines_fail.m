% hydro stream function
clear all
fun=@(x,y,C,a)(C-2*a*y/(x.^2+y.^2-a^2));
a=1;
nC=2;
C=linspace(1,1,nC);
nx=100;
x=linspace(-1,1,nx);
for p=1:nC
    for k=1:nx
        X(p,k)=fzero(@(y)fun(x(k),y,C(p),a),1);
    end
end

for p=1:nC
    plot(x,X(p,:)); hold 'on';
end

%% straight flow + circulation

v0=1;
s2p=2;
dx=0.3;
dy=dx;
x=-2.5:dx:2.5;
y=-2:dy:5;
[X Y]=meshgrid(x,y);

pot=@(a,b,v0,s2p)(-v0*b+s2p*2*log(a.^2+b.^2));
stream=@(a,b,v0,s2p)(-v0*a-s2p*atan(b./a));
Z=pot(X,Y,v0,s2p);
Z2=stream(X,Y,v0,s2p);
[DX,DY]=gradient(Z,dx,dy);
contour(X,Y,Z)
hold on;
contour(X,Y,Z2)
hold on;
quiver(X,Y,DX,DY)
hold off
%% only source

v0=1;
s2p=3;
db=0;
da=2;
dx=0.3;
dy=dx;
edge=5;
x=-edge:dx:edge;
y=-edge:dy:edge;
[X Y]=meshgrid(x,y);

pot=@(a,b,da,db,v0,s2p)(s2p*2*log((a-da).^2+(b-db).^2));
stream=@(a,b,da,db,v0,s2p)(-0*v0*b/10000.5-s2p*atan((b-db)./(a-da)));
Z=pot(X,Y,-da,0,v0,-s2p)+pot(X,Y,da,0,v0,s2p);
Z2=stream(X,Y,-da,db,v0,-s2p)+stream(X,Y,da,db,v0,s2p);
[DX,DY]=gradient(Z);

quiver(X,Y,DX,DY,'linewidth',1);
hold on;
contour(X,Y,Z)
hold on;
contour(X,Y,Z2)
hold off
%%
[x,y] = meshgrid(-2:.2:2, -2:.2:2);
z = x .* exp(-x.^2 - y.^2);
[px,py] = gradient(z,.2,.2);
contour(z), hold on, quiver(px,py), hold off
