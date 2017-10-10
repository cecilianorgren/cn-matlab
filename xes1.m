%% xes1
f=@(v,vt)(exp(-v.^2/vt^2));
df=@(v,vt)(-2*v.*exp(-v.^2/vt^2)/vt^2);
vt=1;
v=0:0.1:2;
ff=f(v,vt);
dff=gradient(ff,v);
dff2=df(v,vt);
plot(v,ff,v,dff,v,dff2)
xlabel('v')
legend('f','df/dv')

%%
wc=0.924;
wh=0.383;
k=wc/vt;
vt=0.9;
wr=0.9*wh;

g2=@(wr,k,vt,wh,wc)(-sqrt(pi/8)*(wh/wc)^2*(wr/k/vt)^3*exp(-0.5*(wr/k/vt)^2));
g2(wr,k,vt,wh,wc)*wr
%(wh/wc)^2
%%
wp=1;
wr=1.2*wp;
k=1;
vt=0.4;
xr=wr/k/vt;
xp=wp/k/vt;
gg=@(xr,xp)(-sqrt(pi)*xp*xr^2*exp(-xr^2));
gg(xr,xp)
%(wh/wc)^2

%%
t=0:10;
g=-1;
f=@(g,t)(exp(2*g*t));
semilogy(t,f(g,t))

%% gamma xes1
t=0:36
t1=2.8; t2=8.1;
U1=0.0041; U2=0.0015;
gamma=log(U2/U1)/2/(t2-t1);
U0=0.0129;
U=@(t,gamma,U0)(U0*exp(2*gamma*t));
plot(t,U(t,gamma,U0))

%% gamma theory
ve=0.4;
k=1;
ope=1;
or=ope*1.2;
gammma=-sqrt(pi)*ope^2*or^2*exp(-(or/k/ve)^2)/(k^3*ve^3)

%%
 n = 2048  nv2 =    0  nlg =    1  mode =    1 
 wp =  1.000   wc =  0.000   qm = -1.000 
 vt1 =  0.000  vt2 =  0.400  v0 =  0.000 
 x1 =  0.100   v1 =  0.000   thetax =  0.000   thetav =  0.000 
 nbins =   50   vlower =  0.000  vupper =  0.000 
%% find max of grwoth rate
gam=@(wr)(wr.^2.*exp(-wr.^2));
wr=0.01:0.01:3;
%plot(wr,gam(wr))
ggam=gam(wr);
grad_gam=gradient(ggam,wr)
plot(wr,gam(wr),wr,grad_gam)
minIndex=find(grad_gam==min(abs(grad_gam)));
[wr(minIndex) grad_gam(minIndex)]
%% find max of grwoth rate
gam=@(x)(sqrt(pi)*(wh/wc/2)^2*x.^3.*exp(-x.^2));
x=0.02:0.01:3;
%plot(wr,gam(wr))
ggam=gam(x);
grad_gam=gradient(ggam,x);
plot(x,gam(x),x,grad_gam); legend('gamma','grad(gamma)')
minIndex=find(grad_gam==min(abs(grad_gam)));
[wr(minIndex) grad_gam(minIndex)]
%%
rw=1;
rr=[0 0.2 0.4 0.6]*rw;
ld=logspace(-0,2,100);
n1=1;
n2=2;
g=1;
x=rw./ld;
do2=@(n1,n2,rw,rr,g,ld)(1+g*2*pi^2*(ld./(rw-rr)).^2.*(n2^2-n1^2));
o2=@(n,rw,rr,g,ld)(1+g*2*pi^2*(ld./(rw-rr)).^2.*(n^2));

fig=figure(1);

for ii=1:4;
    h(ii)=subplot(2,2,ii);
end
for ii=1:4;
semilogx(h(ii),x,do2(1,2,rw,rr(ii),g,ld),...
               x,do2(2,3,rw,rr(ii),g,ld),...
               x,do2(3,4,rw,rr(ii),g,ld),...
               x,do2(4,5,rw,rr(ii),g,ld))
legend(h(ii),'12','23','34','45')
xlabel(h(ii),'r_W/\lambda_D','fontsize',14)
title(h(ii),['\Delta(\omega/\omega_r)^2     (r_R=',num2str(rr(ii)),' r_W)'],...
    'fontsize',14)
end
set(h,'xlim',[x(end) x(1)],'ylim',[0 1e7])
set(gcf,'paperpositionmode','auto')
%%
figure(2);
for ii=1:4;
    h(ii)=subplot(2,2,ii);
end
for ii=1:4;
semilogx(h(ii),x,o2(1,rw,rr(ii),g,ld),...
               x,o2(2,rw,rr(ii),g,ld),...
               x,o2(3,rw,rr(ii),g,ld),...
               x,o2(4,rw,rr(ii),g,ld))
legend(h(ii),'1','2','3','4')
xlabel(h(ii),'r_W/\lambda_D','fontsize',14)
title(h(ii),['(\omega/\omega_r)^2     (r_R=',num2str(rr(ii)),' r_W)'],...
    'fontsize',14)
end
set(h,'xlim',[x(end) x(1)],'ylim',[0 2e7])
set(gcf,'paperpositionmode','auto')

%% bump_on_tail_no.inp
% Parameters
% Mode 1,2,3
k1=1; T1=[1.5 7.9]; vp1=2*pi/diff(T1)/k1;
k2=3; T2=[4.38 8.59]; vp2=2*pi/diff(T2)/k2;
k3=3; T3=[4.22 7.23]; vp3=2*pi/diff(T3)/k3;
%% bump_on_tail_no2.inp
% Parameters
% Mode 1,2,3
k1=1; T1=[0.83 1.44]; vp1=2*pi/diff(T1)/k1;
k2=2; T2=[1.40 2.24]; vp2=2*pi/diff(T2)/k2;
k3=3; T3=[1.13 1.73]; vp3=2*pi/diff(T3)/k3;
[vp1 vp2 vp3]
%% try to sort out relation between number of particles and mass
wp2=0.929;
wp1=0.383;
m1overm2=(wp2/wp1)^2*2^(15-12)
%%
wp2=0.929;
wp1=0.383;
m1=1; m2=m1;
qm1=-1; qm2=qm1;
n1overn2=(wp1/wp2)^2*(qm2/qm1)*(m1/m2)


