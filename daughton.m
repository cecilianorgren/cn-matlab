L=1;
x=-4:0.01:4;
xL=x*L;
irf_units

% Assume B in nT, n in cc, phi in V

% Backgrond magnetic field
Bz0=10;
B=Bz0*tanh(xL);
%
Te=(L*Units.e*Bz0*1e-9).^2*2/Units.mi; % 
Ti=Te;
% Density
n0=Bz0^2*1e-18/8/pi./(Ti+Te); % equilibrium force balance
n=n0.*sech(xL).^2;
%
dB1=0.5*besselj(1,xL*15).*(exp(-abs(xL)));
phiB1=B*1e-9.*dB1*1e-9./(n*1e6)/Units.e/Units.mu0;
%
%f_phiB1=@(deltaB,B,n)(B*1e-9.*deltaB*1e-9./(n*1e6)/Units.e/Units.mu0);

dB2=0.5*besselj(2,xL*15).*exp(-abs(xL));
phiB2=B*1e-9.*dB2*1e-9./(n*1e6)/Units.e/Units.mu0;
%
figure(111)
np=7;
for k=1:np
    h(k)=subplot(np,1,k);
end

plot(h(1),xL,B)
plot(h(2),xL,n)
plot(h(3),xL,Te,'-')
plot(h(4),xL,dB1)
plot(h(5),xL,phiB1)
plot(h(6),xL,dB2)
plot(h(7),xL,phiB2)

%for k=1:np; set(h(k),'xtick',[],'ytick',[]); end 
set(gca,'xtickMode', 'auto');

ylabel(h(1),'B  [nT]')
ylabel(h(2),'n  [cc]')
ylabel(h(3),'Te, Ti')
ylabel(h(4),'dB1')
ylabel(h(5),'phiB1')
ylabel(h(6),'dB2')
ylabel(h(7),'phiB2')

xlabel(h(end),'x/L')

%ylim(h4,[-10 10])