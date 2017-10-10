% goldstein2007
% model
hold off
x=linspace(-1,180,200);
a=14;
ag=10;
phimax=4.8;
dx=107;
phi=@(x,a,phimax,dx)(phimax*sech(x/a-dx/a).^4);
plot(x,phi(x,a,phimax,dx)); hold on;

phigauss=@(x,a,phimax,dx)(phimax*exp(-((x-dx)./a).^2));
plot(x,phigauss(x,ag,phimax,dx),'r'); hold on;

% data
t1=[2007 08 31 10 17 45.50]; t2=[2007 08 31 10 17 45.90]; % one hole
tint=[toepoch(t1) toepoch(t2)];

b=irf_tlim(gsmB4,tint(1),tint(2));
b0=mean(b(:,2:4),1);
b0norm=irf_norm(b0);
%
e3=irf_tlim(gsmE3,tint(1),tint(2));
face3=irf_dot(e3,b0norm);
phi=irf_tappl(irf_integrate(face3-mean(face3(:,2))),'*(1*376)');
plot(phi(:,2),'g')
%b0=