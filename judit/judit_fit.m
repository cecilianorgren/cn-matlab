% Data
Lz = 1;
Lx = [0.5,1,2,5,6,7,8,10,15,20,40];
R = [0.22602171,0.24566273,0.25104229,0.20345172,0.18692082,0.17226913,0.15952006,0.14030094,0.09496913,0.05082819,0.03705537];
% Fit function
xx = linspace(0,40,1000);
Lf = 0.18;
lexp = 0.4;
f0 = 1.8;
f = @(Lx) f0*Lx.*exp(-(Lx/Lf).^lexp);
% Lf = 3;
% lexp = 0.4;
% f0 = 1.5;
% f = @(Lx) f0*(Lx).^(lexp-1)./Lf*lexp.*exp(-(Lx/Lf).^lexp);
Lf = 0.18;
lexp = 0.4;
f0 = 1.6;
exp0 = 0.01;
f = @(Lx) f0*Lx.*(exp(-((Lx/Lf)).^lexp)+exp0./Lx);
L = sqrt(Lx.^2 + Lz.^2);

Lz = 1;
f = @(Lx) Lz./(Lx+3).^0.5-0.15;
f = @(Lx) (Lz./(Lx+5)).^0.5-0.12;
f = @(Lx) (Lz./(Lx)).^0.5-0.15;
%f = @(Lx) atan((Lz)./(Lx)+0.02).^0.95;
%f = @(Lx) atan(2*(Lz)./(Lx)+1.00).^1.0;

%Lz/Lx

plot(Lx,R,'*',xx,f(xx),'-')
%plot(L,R,'*')
%semilogx(Lx,R,'*',xx,f(xx),'-')
hca = gca;
hca.YLim = [0 0.4];
grid on;

%%
% Judits fit to data 
eta = logspace(-2,1,1000);
a = 0.18435;
b = 0.30802;
c = 0.02164;
fR = @(eta) a*eta./(b+eta)+c;
plot(1./eta,fR(eta),'*')