clear all
f1=@(z,p,L)(((z-L)./p)./sqrt(p.^2+(z-L).^2));
L=10;
nz=100;
np=100;
z=linspace(0,L,nz);
p=logspace(-4,log10(5*L),np);
%p=linspace(0.01,50,np);
%plot(z,f1(z,p,L)-f1(z,p,0))


for a=1:np
    f(a,:)=f1(z,p(a),L);
end
ax=pcolor(z,p,f)
ac=colorbar;
caxis([-8 0.25])
shading('flat')
ylabel('rho')
xlabel('z')

%% integral x^2/(x^2+a^2)
a=10;
fun=@(x,a)(x.^2./(x.^2+a.^2));

nx=1000;
r=linspace(1e-3,1e4,nx);

na=100;
a=linspace(2e-3,1e1,na);
clear a
a=10;

for k=10:nx
    y=fun(r(1:k),a);
    int(k-9)=trapz(r(1:k),y);
end

answer=@(x,a)(-a.*tan(x./a).^(-1));
%%
fun=@(x,a)(x.^2./(x.^2+a^2));
na=100;
a=linspace(0.01,10,na);

%%
nr=1000;
rmax=200;
r=linspace(0.0,rmax,nr);
rpmax=10;
rp=0:0.001:rpmax;
answer=@(x,a)(a-x.*atan(a./x));
for k=1:nr;
answ(k)=trapz(rp,fun(rp,r(k)));
answan(k)=answer(r(k),rpmax);
end
plot(r,answ,rpmax,max(answ),'*')
xlabel('r');ylabel('pot')

%%
% draw potentials from charges
irf_units;
phi=@(r,qs)(qs*Units.e./(4*pi*Units.eps0*r.^2));
phi=@(r,qs)(qs./(r.^2));
if 0
r=linspace(1e-2,50,100);
theta=linspace(0,2*pi,100);
x=r.*cos(theta);
y=r.*sin(theta);
else
    x=logspace(-3,-1.5,50);
    x=[-x(end:-1:1) x];
    y=x;
end

[X,Y]=meshgrid(x,y);


%R = sqrt(X.^2 + Y.^2) + eps;
%Z = sin(R)./R;
x1=0;
x2=1e-2;
Z=phi(sqrt((X-x1).^2+Y.^2),-1)+phi(sqrt((X-x2).^2+Y.^2),-1);
[DX,DY]=gradient(Z,1e-1,1e-1);
contour(X,Y,Z);
quiver(X,Y,DX,DY);
shading flat
%% example mesh
figure;
[X,Y] = meshgrid(-8:.5:8);
R = sqrt(X.^2 + Y.^2) + eps;
Z = sin(R)./R;
mesh(Z);